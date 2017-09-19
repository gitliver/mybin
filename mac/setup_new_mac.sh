#!/bin/bash

# A script for installing programs on a new Mac
# June 1, 2017

echo "*** install brew ***"
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

for i in wget coreutils tmux youtube-dl imagemagick gpg python3 postgres; do
	echo "*** install "${i}" ***"; 
	brew install ${i}; 
done

echo "*** install pip ***"
wget https://bootstrap.pypa.io/get-pip.py
sudo python get-pip.py

echo "*** install virtualenv ***"
sudo pip install virtualenv
