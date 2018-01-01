#!/bin/bash
ffmpeg -f image2 -r 8 -i image%05d.png -vcodec mpeg4 -y thermal.mp4
