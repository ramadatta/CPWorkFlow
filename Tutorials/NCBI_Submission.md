## Raw Sequence Reads and Genome Assembly Submission

* Log into NCBI account with credentials
* We need to submit the data to Biosample first and then use Biosample's Accession to create a Bioproject Accession

### BioSample Submission
* While submitting BioSample, we need to fill Microbe1.0 template which required the following fields required to submit (Green are madatory): 
![image](https://user-images.githubusercontent.com/3212461/125230295-1d42e200-e30b-11eb-8d30-bb3e73b1ef88.png)![image](https://user-images.githubusercontent.com/3212461/125230371-39468380-e30b-11eb-9bae-a794c36ce2d9.png)![image](https://user-images.githubusercontent.com/3212461/125230386-43688200-e30b-11eb-97e1-a54d2102ac79.png)

And any one of the blue is required ![image](https://user-images.githubusercontent.com/3212461/125230457-609d5080-e30b-11eb-9d97-9cd447727327.png)

### Solve the "multiple BioSamples cannot have identical attributes" error when uploading data to NCBI!
_ERROR:
Your table upload failed because multiple BioSamples cannot have identical attributes. You should have one BioSample for each specimen, and each of your BioSamples must have differentiating information (excluding sample name, title, bioproject accession and description). This check was implemented to encourage submitters to include distinguishing information in their samples. If the distinguishing information is in the sample name, title or description, please recode it into an appropriate attribute, either one of the predefined attributes or a custom attribute you define. If it is necessary to represent true biological replicates as separate BioSamples, you might add an 'aliquot' or 'replicate' attribute, e.g., 'replicate = biological replicate 1', as appropriate. Note that multiple assay types, e.g., RNA-seq and ChIP-seq data may reference the same BioSample if appropriate._

# Answer1
What should I do if I encounter this problem when uploading 16s/RNAseq/chip_seq?

After searching all night, I found that almost no one was talking about the point!

The key question: In the
table, it 绿色is 必填项necessary to ensure that at least one factor can distinguish each sample (except the name) . In other words, in addition to the sample name, one of the remaining mandatory items must be different in each sample.
Solution:
If you, like me, except that the sample is not the same as the name of the other information is not, then I recommend you directly modify *lat_lon, which is the latitude and longitude . Since the latitude and longitude are kept to two decimal places, we directly add 0.01 at a time ( just drop down in EXCEL to fill it quickly! ), which has almost no effect on the position, and the samples can also be distinguished, which is a perfect solution.
If you have any questions, directly send a private message to the background.

# Answer2
In my case, I just pasted all the *sample_name* into the *isolate* column. This solved the problem!
