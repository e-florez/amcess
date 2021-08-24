cd ../p_ascec
gfortran -o ascec-v03 main_ascec.f config.f function_E.f ran0.f draw_molekel.f masscen.f symelem.f wtelem.f
cp ascec-v03 ../
mv ascec-v03 ../example
