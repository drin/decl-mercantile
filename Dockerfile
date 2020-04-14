FROM ubuntu:18.04

# Popper Ceph ARGs
ARG DEBIAN_FRONTEND=noninteractive
ARG GIT_URL="https://github.com/uccross/skyhookdm-ceph.git"
ARG GIT_REF="skyhook-luminous"
ARG PATH_TO_REPO="/code/skyhookdm-ceph"
ARG EXTRA_PKGS=""

# ------------------------------
# Grab skyhook source and prep system dependencies
RUN    apt-get update                                                             \
    && apt-get install -y                                                         \
               git                                                                \
               gnupg2                                                             \
               ccache                                                             \
               sudo                                                               \
               software-properties-common                                         \
               dirmngr                                                            \
               apt-transport-https                                                \
               lsb-release                                                        \
               ca-certificates                                                    \
    && sh -c 'if [ -n "$EXTRA_PKGS" ]; then apt-get install -y "$EXTRA_PKGS"; fi' \
    && git clone --branch $GIT_REF --depth 1 $GIT_URL $PATH_TO_REPO               \
    && cd $PATH_TO_REPO                                                           \
    && ./install-deps.sh

# ------------------------------
# Compile simple skyhook targets
# NOTE: this is basically in `entrypoint.sh`
# RUN    cd $PATH_TO_REPO     \
#     && ./do_cmake.sh        \
#     && cd build             \
#     && make -j2 cls_tabular \
#     && make -j4 run-query   \
#     && echo `date`

# Personal ARGs
ARG CONFIGS_GIT_URL="https://github.com/drin/configs.git"
ARG PATH_TO_CONFIGS="/code/configs"

# ------------------------------
# Prep personal environment
RUN    git clone --branch "cloudlab" $CONFIGS_GIT_URL $PATH_TO_CONFIGS   \
    && cd $PATH_TO_CONFIGS                                               \
    && CONFIG_ROOT=$(pwd) bash installers/setup-environment.bash         \
    && CONFIG_ROOT=$(pwd) bash installers/install-tools.bash             \
    && CONFIG_ROOT=$(pwd) bash installers/install-alacritty.ubuntu.bash  \
    && CONFIG_ROOT=$(pwd) bash installers/install-vim-plugins.bash setup \
    && apt-get clean -y                                                  \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* debian/

ARG CODE_GIT_URL="https://github.com/drin/decl-mercantile.git"
ARG PATH_TO_CODE="/code/decl-mercantile"
RUN git clone --branch "cloudlab" $CODE_GIT_URL $PATH_TO_CODE

COPY entrypoint.sh /

ENTRYPOINT ["/entrypoint.sh"]

