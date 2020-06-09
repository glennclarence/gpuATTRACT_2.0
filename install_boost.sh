sudo apt-get update
#Install required packages
sudo apt-get -y install build-essential g++ python-dev autotools-dev libicu-dev libbz2-dev

mkdir boost
wget http://downloads.sourceforge.net/project/boost/boost/1.54.0/boost_1_54_0.tar.gz
tar -zxvf boost_1_54_0.tar.gz
rm -f boost_1_54_0.tar.gz
cd boost_1_54_0

#Compile and install
./bootstrap.sh --prefix=/boost
cpuCores=`cat /proc/cpuinfo | grep "cpu cores" | uniq | awk '{print $NF}'`
echo "Available CPU cores: "$cpuCores
sudo ./b2 --with=all -j $cpuCores install
cd ..
rm -rf boost_1_54_0
