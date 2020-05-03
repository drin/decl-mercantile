"""
Module that contains various connectors. The most important, for now, is a rados connector.

This module is necessary to provide an insulating layer above dependencies, so that this library
can (hopefully) be agile.
"""

import os
import sys
import logging

# functions
from skyhookdm_singlecell.util import try_import

# TODO: this is temporary to see if the ubuntu package supports a decent version of librados
sys.path.insert(0, '/usr/lib/python3/dist-packages')

# optional, dynamic module imports
rados = try_import('rados')


class RadosIOContext(object):
    @classmethod
    def for_cluster_pool(cls, cluster_conn, ceph_pool):
        return cls(cluster_conn.open_ioctx(ceph_pool))

    def __init__(self, rados_io_ctx, **kwargs):
        super().__init__(**kwargs)

        self.io_context = rados_io_ctx

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.io_context.close()

    def write_data(self, storage_obj_name, storage_obj_data):
        self.io_context.aio_write_full(storage_obj_name, storage_obj_data)
        self.io_context.set_xattr(
            storage_obj_name,
            'size',
            str(len(storage_obj_data)).encode('utf-8')
        )


class RadosConnector(object):
    """
    Class that tries to provide a usable, composable interface over the rados library.

    For details on the rados library, see:
        https://docs.ceph.com/docs/master/rados/api/python/

    For details on how to install the rados library, see ubuntu 18.04 package: 
        python3-rados_12.2.12-0ubuntu0.18.04.5_amd64.deb

    """

    logger = logging.getLogger('{}.{}'.format(__module__, __name__))
    logger.setLevel(logging.INFO)

    # explicit configuration, currently tailored to skyhook cloudlab installation script defaults
    default_config = os.path.join(os.environ['HOME'], 'cluster', 'ceph.conf')
    default_conf_overrides = {
        'keyring': os.path.join(os.environ['HOME'], 'cluster', 'ceph.client.admin.keyring')
    }

    @classmethod
    def check_for_rados(cls):
        """
        Class method for convenience that checks if rados was successfully imported.
        """

        if rados is None:
            err_msg = 'Could not import python module: "rados"'

            cls.logger.error(err_msg)
            sys.exit(err_msg)

    @classmethod
    def connection_for_config(cls, path_to_config=default_config, **kwargs):
        """
        Class method to configure a cluster handle and initialize a RadosConnector instance with
        the configured handle.

        For official Ceph documentation, see:
            https://docs.ceph.com/docs/master/rados/api/python/#configure-a-cluster-handle
        """

        cls.check_for_rados()

        # make local copy of default overrides
        conf_overrides = dict(cls.default_conf_overrides)

        # update local copy of overrides based on passed keyword arguments
        for kw_name, kw_val in kwargs.items():
            if kw_name in conf_overrides:
                conf_overrides[kw_name] = kw_val

        # create a cluster connection using a config file and configuration overrides
        cluster_conn = rados.Rados(conffile=path_to_config, conf=conf_overrides)

        cls.logger.info(f'librados version: {str(cluster_conn.version())}')

        return cls(cluster_conn)

    def __init__(self, cluster_conn, **kwargs):
        super().__init__(**kwargs)

        self.cluster_conn = cluster_conn
        self.is_connected = False

    def cluster_info(self):
        if not self.is_connected:
            return 'Not yet connected to cluster'

        return 'Cluster ({}) statistics: {}'.format(
            self.cluster_conn.get_fsid(),
            '\n'.join([
                f'\t{stat_key:12s} => {stat_val}'
                for stat_key, stat_val in self.cluster_conn.get_cluster_stats().items()
            ])
        )

    def connect(self):
        self.__class__.logger.info('Attempting connection to: {}'.format(
            self.cluster_conn.conf_get('mon initial members')
        ))

        self.cluster_conn.connect()
        self.is_connected = True

        return self

    def context_for_pool(self, pool_name):
        return RadosIOContext.for_cluster_pool(self.cluster_conn, pool_name)


if __name__ == '__main__':
    print('Success')
