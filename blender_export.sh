#!/bin/sh
if [ -z "$2" ] ; then
	echo "usage: $0 <input.blend> <output.ecs>" > /dev/stderr
	exit 1
fi
blender -noaudio -b "$1" -P "$(dirname $0)/blender_export.py" -- "$2"
