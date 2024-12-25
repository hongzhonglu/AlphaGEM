#./bash

tar -xyz data_source.tar.xz 


cd ./tools/eggnog_mapper
python download_eggnog_data.py
echo y
echo y
echo y

cd ./tools/clean
python build.py install

cd ./tools/OrthoFinder
python setup.py




