#!/bin/bash
cd ~/singularity-ce-4.0.1

#Needed to mount google drive before running this

#sudo mkdir /mnt/g
#sudo mount -t drvfs G: /mnt/g  

export SINGULARITY_BIND='/mnt/g/.shortcut-targets-by-id/0B3ejcGO9fA5eYXlHMW9tekVERmc/JB labs documents/Updated_JB_documents_2019/JB LAB_Projects/PUG_PPIs':/home/heineike/PUG_PPIs,'/mnt/c/Users/heine/OneDrive/Documents/GitHub':/home/heineike/github
