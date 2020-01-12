# DeclMercantile
A repository for declarative programmable storage work by Aldrin. I named this repository
decl-mercantile as a related, but different, project with respect to [DeclStore][paper-declstore].

# Overview
There has been some flux in the task priority for declarative programmable storage. This project
(or at least, the part that this repository maintains) will first focus on how existing gene
expression data can be serialized in a format for skyhook to ingest, and then the performance
profiles of how skyhook can query the data.

# Architecture
The directory layout of this repository is as follows (at the time of this writing):

    3 directories, 3 files
    .
    ├── code
    │   ├── poetry.lock
    │   ├── pyproject.toml
    │   ├── README.rst
    │   ├── skyhookdm_singlecell
    │   ├── tests
    │   └── toolbox
    ├── LICENSE
    ├── README.md
    ├── scripts
    │   └── shell-utils
    ├── setup.sh
    └── submodules
        └── ceph-container

    8 directories, 6 files

### Code

The _code_ directory contains a [poetry][tool-poetry] project, ``skyhookdm_singlecell``. This
project currently contains simple code for:

* Parsing some [HCA data][project-hca]
* Parsing another simple gene expression format
* Serializing gene expression data in [Arrow][project-arrow] format with extra
  [Skyhook][project-skyhook] metadata
    * This code is python and leans heavily on [pyarrow][lib-pyarrow] for serializing in Arrow
      format.

The code will eventually expand to include anything lower-level related to my research and
declarative programmable storage. For anything higher-level, the repository [XHCA][lib-xhca] will
be used.

### Submodules

Though a work in progress (especially since I have little experience with submodules), this is
where git submodules are. Initially, I wanted to have ceph be a submodule so that this repository
could house any modifications to ceph and still be buildable, but I have not pushed that aspect
forward in awhile.

### Scripts

Here you may find standalone scripts that provide some independent functionality. At the moment,
there is a script based on [xweichu's][developer-aaron] ceph build script for cloudlab. I have not
tested my version (scripts/shell-utils/run-ceph-container.fish), but that will hopefully happen in
the near-ish future.


<!-- Resources -->
[paper-declstore]: https://www.usenix.org/conference/hotstorage17/program/presentation/watkins

[tool-poetry]:     https://python-poetry.org/

[project-arrow]:   https://arrow.apache.org/
[project-skyhook]: https://sites.google.com/view/skyhook-programmable-storage
[project-hca]:     https://data.humancellatlas.org/

[lib-pyarrow]:     https://arrow.apache.org/docs/python/
[lib-xhca]:        https://github.com/disorderlylabs/xhca

[developer-aaron]: https://github.com/xweichu
