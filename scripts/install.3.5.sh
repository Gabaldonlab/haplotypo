#!/bin/bash
###
# Haplotypo installer for UNIX.
# version 0.1a
###

branch="master" 
if [ ! -z $1 ]; then branch=$1; fi

echo "##################################################################################################"
echo "#                                                                                                #"
echo "#                                      Haplotypo installer                                       #"
echo "#                                                                                                #"
echo "#        version 0.1a                                                                            #"
echo "##################################################################################################"
echo ""
echo "Haplotypo and its dependencies will be installed in:" `pwd`/haplotypo
echo "Installation will take 5-10 minutes. "
echo ""

# sleep
echo "I'll proceed with installation in 2 seconds... Press Ctrl-C to cancel."
sleep 2s

echo ""
echo `date` "Checking dependencies..."

exists()
{
  command -v "$1" >/dev/null 2>&1
}

error=""
# check if all basic programs exists
for cmd in echo awk git wget unzip tar nano gcc g++ make cd ln date ldconfig pip; do
    if ! exists $cmd; then
        case $cmd in
            "pip")
                echo "Install pip first (ie. 'sudo apt-get install python-pip')!"  
                ;;
            *)
                echo "Install $cmd first (ie. 'sudo apt-get install $cmd')!"
                ;;
        esac
        error=1
    fi
done

# check if all libs present #BWA
for lib in libz; do
    if [ -z "$(ldconfig -p | grep $lib.so)" ] ; then
        echo " Missing library $lib !"
        error=1
    fi
done

# check headers #BWA
for lib in zlib.h; do
    if [ ! -s /usr/include/$lib ] && [ ! -s /usr/lib/$lib ]; then
        echo " Missing headers $lib !"
        error=1
    fi
done

# skip if error
if [ ! -z $error ]; then
    echo -n "We've finded missing dependencies.\n"
    echo -n "Would you like to install (y or n)?\n"
    read response
    if [ $response != "y" ]; then        
        echo "\nAborted due to missing dependencies (see above)!"
        return 1;
    else
        sh install_dependencies.sh
    fi
fi

# check python version 
PyVer=`python --version 2>&1 | cut -f2 -d" " | cut -f-2 -d"."`
if [ $PyVer != "2.7" ] && [ $PyVer != "2.6" ]; then 
    echo ""
    echo "[ERROR] Install Python 2.7!"
    echo "If you have Python 2.7 already installed, you can either "
    echo "make an alias before installation and use of Redundans ('alias python=python2.7' should do)"
    echo "or use Python virtual environment (https://virtualenv.pypa.io)."
    return 1
fi
echo " Everything looks good :) Let's proceed..."
sleep 2s

pwd=`pwd`

'''
Dependencies versions
'''
GATK_VERSION=4.0.12.0
HTSLIB_VERSION=1.9
SAMTOOLS_VERSION=1.9
BCFTOOLS_VERSION=1.9
BWA_VERSION=0.7.15

echo " Creating dependencies folder..."
mkdir dependencies
cd dependencies

echo " Installing basic software..."
apt-get install -y software-properties-common
apt-get install -y build-essential

echo " Installing Java..."
add-apt-repository -y ppa:webupd8team/java
apt-get update
echo debconf shared/accepted-oracle-license-v1-1 select true | debconf-set-selections
echo debconf shared/accepted-oracle-license-v1-1 seen true | debconf-set-selections
apt-get install -y --force-yes oracle-java8-installer

echo "Installing BWA"
wget https://github.com/lh3/bwa/releases/download/v$BWA_VERSION/bwa-$BWA_VERSION.tar.bz2
tar -vxjf bwa-$BWA_VERSION.tar.bz2
cd bwa-$BWA_VERSION
make 
cd ..

echo "Installing Samtools, Bcftools and Htslib..."
apt-get update
apt-get install -y libbz2-dev
apt-get install -y bzip2
apt-get install -y zlib1g-dev
apt-get install -y libncurses5-dev 
apt-get install -y libncursesw5-dev
apt-get install -y liblzma-dev

wget https://github.com/samtools/htslib/releases/download/$HTSLIB_VERSION/htslib-$HTSLIB_VERSION.tar.bz2
tar -vxjf htslib-$HTSLIB_VERSION.tar.bz2
cd htslib-$HTSLIB_VERSION
make
cd ..

wget https://github.com/samtools/samtools/releases/download/$SAMTOOLS_VERSION/samtools-$SAMTOOLS_VERSION.tar.bz2
tar -vxjf samtools-$SAMTOOLS_VERSION.tar.bz2
cd samtools-$SAMTOOLS_VERSION
make
cd ..

wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-$BCFTOOLS_VERSION.tar.bz2
tar -vxjf bcftools-$BCFTOOLS_VERSION.tar.bz2
cd bcftools-$BCFTOOLS_VERSION
make
cd ..

echo "Installing GATK..."
wget https://github.com/broadinstitute/gatk/releases/download/$GATK_VERSION/gatk-$GATK_VERSION.zip
unzip gatk-$GATK_VERSION.zip

echo "Installing Python packages"
apt-get -y install python3-pip
pip3 install numpy
pip3 install biopython
pip3 install psutil
pip3 install pysam
pip3 install pyvcf
python -m pip3 install --user matplotlib ipython jupyter pandas sympy nose seaborn


dep_folder=/home/dependencies
PATH="$dep_folder/samtools-1.9/:${PATH}"
PATH="$dep_folder/bcftools-1.9/:${PATH}"
PATH="$dep_folder/bwa-0.7.15/:${PATH}"
echo 'alias haplotypo="python $(pwd)/bin/3.5/haplotypo.py"' >> ~/.bashrc

apt-get clean
set -x; rm -rf /var/lib/apt/lists/*

echo `date` "Installation finished!"
echo "##################################################################################################"
echo "# Haplotypo depends on several programs (https://github.com/Gabaldonlab/haplotypo#prerequisities)#"
echo "# Acknowledge their authors, get familiar with licensing and register if necessary.              #"
echo "##################################################################################################"
echo ""
