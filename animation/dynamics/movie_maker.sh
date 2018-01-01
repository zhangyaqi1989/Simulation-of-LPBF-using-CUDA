#!/bin/bash
ffmpeg -f image2 -r 8 -i image%03d.png -vcodec mpeg4 -y dynamic.mp4
