#!/usr/bin/env bash
DEEP_SCAN=${DEEP_SCAN:=false}

echo "========================================"
echo "Clearing runner disk space"
echo "========================================"

echo "Listing 20 largest packages"
dpkg-query -Wf '${Installed-Size}\t${Package}\n' | sort -n | tail -n 20

if DEEP_SCAN; then
  echo
  echo "Running deep scan"
  du -h / 2>/dev/null | grep "^ *[0-9]*G"
fi

echo
echo "Disk space usage: before"
df -h

sudo rm -rf /opt/az
sudo rm -rf /usr/lib/google-cloud-sdk
sudo rm -rf /usr/lib/jvm
sudo rm -rf /opt/google/chrome
sudo rm -rf /usr/lib/firefox
sudo rm -rf /opt/microsoft/powershell
sudo rm -rf /usr/share/dotnet
sudo rm -rf /opt/ghc
sudo rm -rf /opt/hostedtoolcache

echo
echo "Disk space usage: after"
df -h
