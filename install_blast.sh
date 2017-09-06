#!/usr/bin/env bash

#ead  -n 1 -p "Input Selection:" mainmenuinput
#tar zxvpf ncbi-blast+2.2.29-x64-linux.tar.gz
#open "ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/"

FILE=$1
TARGETDIR=$2
DIRNAME="$(dirname $FILE)"
echo "************ INSTALLING BLAST ************"
echo "File:             $FILE"
echo "File Dir Path:    $DIRNAME"
echo "Curr Dir:         $PWD"
echo "Target Dir:       $TARGETDIR"
echo "******************************************"
cd $DIRNAME
if [ ! -f $FILE ]; then
    echo "File $FILE not found!"
    echo "Run open "ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/""
    echo "Login as guest."
    echo "Download blast package as a tar"
    exit 1
fi



# Extract
o="blast_temp_install"
rm -r $o
mkdir $o
echo "  ...extracting blast"
tar zxvpf $FILE -C $o &>/dev/null
echo "     > blast extracted"
EXTRACTED="$(ls $o)"

# Move to target dir
DEST="$TARGETDIR/$EXTRACTED"
echo "BLAST may already b installed at $DEST"
echo -n "Cancel and overwrite current installation? (y/n): "
if [ -d "$DEST" ]; then
    read answer
    if echo "$answer" | grep -iq "^y" ;then
        echo "User continued installation."
        rm -r $DEST
    else
        echo "Installation canceled."
        exit 1
    fi
fi
mv "$o/$EXTRACTED" "$DEST"
echo "   ...Moving extracted $mv"

# Export path
newpath=targetdir/$EXTRACTED/bin
export PATH=$PATH:$newpath
echo "    ...Saving path to bash profile: $TARGETDIR/$newpath"
echo $PATH

# Create new dir
NEWDIR="$TARGETDIR/blastdb"
echo "   ...Making new folder $NEWDIR"
echo "$NEWDIR already exists."
echo -n "Overwrite $NEWDIR? (y/n): "
if [ -d "$NEWDIR" ]; then
    read answer
    if echo "$answer" | grep -iq "^y" ;then
        echo "User continued installation."
        rm -r $NEWDIR
    fi
fi

mkdir $NEWDIR
export BLASTDB=$NEWDIR
cd $TARGETDIR
#ls
echo "******* BLAST installation complete *******"