# pca_geophysical_fields
Code available to calculate the principal components between geophysical fields. 

The code is focused in the Chilean subduction zone. It is written in Jupyter interface using python 2.7
and you just need to be concern about two things for reproduce the plots and results shown in the paper.

                                               
The main code that you need to use is main_pca and it is formed by 15 cells. You need to use the libraries
shown in cell two in order to extract the principal components.

First. Make sure that the folders  "data_fields" and "data grids" are in the same path of the scripts.
These two folders contain the file inputs (gravity, friction and locking grids. ) to calculate the principal components.

Second. You will need to give to the code some input parameters along the procedure. 
1th input in cell 2: The path of work. Is the directory where you have the scripts and data -> example = /user/desktop/pca/scripts_data/
/script_data is where you have the files
2th input in cell 4: The number of perpendicular profiles and  then length of them in meter units. example-> 144 ... 150000
                     Therefore you want to use 144 profiles of 150 km each one
3th input in cell 8: Here 3 inputs are needed. The number of fields to analize (c). example-> 3
                                               The name of field to analize but with number nomenclature. For gravity = 1, locking = 2 and friction = 3
                                                         So if you want to analize the three fields you need to enter 1,2,3
                                                         If you want to analize gravity and friction you will need enter 1,3
                                               The mode you want to study (k), if you want to study the first mode , enter 1
                                               
note : It might be that you need to install the package -forge basemap data hire- even if you have basemap already installed.
to install run in the terminal --> conda install -c conda-forge basemap-data-hires
                                            

