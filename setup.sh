#/bash



cd ./tools/eggnog_mapper
python download_eggnog_database.py
echo y
echo y
echo y

cd ./tools/clean
python build.py install


