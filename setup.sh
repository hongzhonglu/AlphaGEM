#/bash



cd ./tools/eggnog_mapper
python download_eggnog_data.py
echo y
echo y
echo y

cd ./tools/clean
python build.py install


