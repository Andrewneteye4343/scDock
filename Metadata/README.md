If you want to group samples during input, please provide a metadata.txt file following the demonstrated format.  
When metadata.txt is provided, the Seurat object will include a metadata column named "sample_group", defined according to the file.  
The metadata.txt file should contain two tab-separated columns:  
● file_name  
● group
For example, with two input files (sample1.txt and sample2.txt):  
In this case, sample1 belongs to the Control group, while sample2 belongs to the Treatment group.
