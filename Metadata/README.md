If you want to group samples during input, please provide a metadata.txt file following the demonstrated format.  
When metadata.txt is provided, the Seurat object will include a metadata column named "sample_group", defined according to the file.  
The "sample_group" label can be used for argument "Run_CellChat_group_by" in CellChat analysis.  
The metadata.txt file should contain two tab-separated columns:  
● file_name  
● group  

```text
file_name    group
sample1.txt  Control
sample2.txt  Treatment
```

For example, with two input scRNA files (sample1.txt and sample2.txt):  
In this case, sample1 belongs to the Control group, while sample2 belongs to the Treatment group.
