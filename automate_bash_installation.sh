#!/usr/bin/env bash

hostname="ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST"

ftp -d -in $hostname << ftpEOF
quote USER "anonymous"
ftpEOF