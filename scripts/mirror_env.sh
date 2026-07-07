#!/usr/bin/env bash
# Shared connection + runtime settings for the ChromBERT server mirror.
# Source this from the other scripts:  source "$(dirname "$0")/mirror_env.sh"
#
# The mirror does NOT persist changes unless the user snapshots it, and the login
# PASSWORD rotates on restart — so we authenticate with an SSH KEY (persisted in
# the mirror's /root/.ssh/authorized_keys, which lives under the persistent /root).
# The IP is FIXED; only the PORT may change on restart. Override MIRROR_PORT per
# session if it moved:  MIRROR_PORT=12345 ./scripts/run_chrombert_region_emb.sh ...

: "${MIRROR_IP:=172.16.78.10}"     # fixed
: "${MIRROR_PORT:=35963}"          # may change on restart — override via env
: "${MIRROR_USER:=root}"
: "${MIRROR_KEY:=$HOME/.ssh/ecr_navigator}"

# Host-key pinning is disabled: the mirrors are ephemeral and reuse the same IP
# with different host keys per instance, so a pinned known_hosts entry just causes
# "Host key verification failed". These are trusted internal hosts on a private IP.
HOSTKEY_OPTS="-o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null"
SSH="ssh -i $MIRROR_KEY $HOSTKEY_OPTS -p $MIRROR_PORT ${MIRROR_USER}@${MIRROR_IP}"
SCP_BASE="scp -i $MIRROR_KEY $HOSTKEY_OPTS -P $MIRROR_PORT"   # note: scp uses -P (uppercase)

# Run a command inside the ChromBERT conda env with the required env vars.
# ChromBERT import needs PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION=python (stale
# protobuf/onnx). HF_ENDPOINT points at a reachable mirror for downloads.
mirror_py () {
  $SSH "source /opt/conda/etc/profile.d/conda.sh; conda activate base; \
        export PATH=/root/bin:\$PATH \
               PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION=python \
               HF_ENDPOINT=${HF_ENDPOINT:-https://hf-mirror.com}; \
        $*"
}
