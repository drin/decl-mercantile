function find-docker-network
    if test -z $argv[1]
        echo "Usage: find-docker-network <network-name>"
        return 1
    end

    docker network ls -f name=$argv[1] | tail -1 | xargs echo | cut -d ' ' -f 2
end

function create-ceph-network
    docker network create --subnet=113.1.1.240/28 ceph-network
end

function debug-ceph-mds
    docker run -it --name=mds                                                            \
                   --net=ceph-network                                                    \
                   --ip 113.1.1.244                                                      \
                   -v /raiddata/research-mounts/by-server/mon/etc/ceph:/etc/ceph         \
                   -v /raiddata/research-mounts/by-server/mon/var/lib/ceph:/var/lib/ceph \
                   -v /raiddata/research-mounts/by-server/mon/var/log/ceph:/var/log/ceph \
                   -e CEPHS_CREAT=1                                                      \
                   ceph/daemon mds
end

function runc-ceph-mon
    if ! test (find-docker-network ceph-network) = "ceph-network"
        echo "Please create the virtual docker network 'ceph-network' first"
        return 1
    end

    docker run -d --name=mon                                                            \
                  --net=ceph-network                                                    \
                  --ip 113.1.1.242                                                      \
                  -v /raiddata/research-mounts/by-server/mon/etc/ceph:/etc/ceph         \
                  -v /raiddata/research-mounts/by-server/mon/var/lib/ceph:/var/lib/ceph \
                  -v /raiddata/research-mounts/by-server/mon/var/log/ceph:/var/log/ceph \
                  -e MON_IP=113.1.1.242                                                 \
                  -e CEPH_PUBLIC_NETWORK=113.1.1.240/28                                 \
                  ceph/daemon mon

    docker run -d --name=mgr                                                            \
                  --net=ceph-network                                                    \
                  --ip 113.1.1.243                                                      \
                  -v /raiddata/research-mounts/by-server/mon/etc/ceph:/etc/ceph         \
                  -v /raiddata/research-mounts/by-server/mon/var/lib/ceph:/var/lib/ceph \
                  -v /raiddata/research-mounts/by-server/mon/var/log/ceph:/var/log/ceph \
                  -e CEPHS_CREAT=1                                                      \
                  ceph/daemon mgr

    docker run -d --name=mds                                                            \
                  --net=ceph-network                                                    \
                  --ip 113.1.1.244                                                      \
                  -v /raiddata/research-mounts/by-server/mon/etc/ceph:/etc/ceph         \
                  -v /raiddata/research-mounts/by-server/mon/var/lib/ceph:/var/lib/ceph \
                  -v /raiddata/research-mounts/by-server/mon/var/log/ceph:/var/log/ceph \
                  -e CEPHS_CREAT=1                                                      \
                  ceph/daemon mds
end

function runc-ceph-osd
    if test -z $argv[1]; or test -z $argv[2]
        echo "Usage: runc-ceph-osd <osd-name> <osd-device>"
        return 1
    end

    set osd_name     $argv[1]
    set osd_dev_path $argv[2]

    set osd_bindmount "/raiddata/research-mounts/by-server/$osd_name"

    if ! test -d $osd_bindmount
        echo "Mounts for OSD '$osd_name' not found."
        echo "Check dir '$osd_bindmount'"
        return 1
    end

    if ! test -b $osd_dev_path
        echo "Device '$osd_dev_path' is invalid."
        return 1
    end

    docker run -d --name=$osd_name                             \
                  --net=ceph-network                           \
                  --device=$osd_dev_path                       \
                  -v $osd_bindmount/etc/ceph:/etc/ceph         \
                  -v $osd_bindmount/var/lib/ceph:/var/lib/ceph \
                  -v $osd_bindmount/var/log/ceph:/var/log/ceph \
                  -v $osd_dev_path:/dev/sdh                    \
                  -e OSD_DEVICE=/dev/sdh                       \
                  ceph/daemon osd_ceph_disk
end
