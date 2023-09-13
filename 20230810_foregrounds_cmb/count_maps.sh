for f in output/galactic/*
do
    echo $f, $(ls $f | wc -l)
done
