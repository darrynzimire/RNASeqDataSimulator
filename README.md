# RSDS: A simulator for raw RNA-sequencing data 

The RNA-seq data simulator (RSDS) is a command-line interface implemented in Python 3. The tool simulates raw RNA-sequencing data by emulating characteristics of real RNA-seq data. Parameters to control the properties of the simulated data are available as tuneable settings such as fragment-length distribution, customized Phred-quality score modelling, customized transcript expression profiling, as well as differential transcript expression simulation. 

### Prerequisites

* Python 3
* Pyfaidx
* Numpy

## Usage

Here is a list of all the parameters for running the RSDS program

| Parameters  | Description                                                                       |
|-------------|-----------------------------------------------------------------------------------|
| -r          |  Read length (default value = 100)                                                |
| -n          | Number of reads to simulate                                                       |
| -f          | Reference transcriptome file in FASTA format                                      |
| -s          | Random seed value for reproducibility                                             |
| -o          | Output file prefix                                                                |
| -seqmodel   | Phred quality-score model                                                         |
| -countModel | Read count table model file                                                       |
| -FL         | Fragment length distribution parameters (default: mean=250, standard deviation=25 |
| -SE         | Single-end RNA-seq data                                                           |
| -PE         | Paired-end RNA-seq data                                                           |


## Phred quality-score model

A list of the parameters for generating a quality-score model 

| Parameters | Description                        |
|------------|------------------------------------|
| -i         | FASTQ file (read 1)                |
| -i2        | FASTQ file (read 2)                |
| -o         | Output file prefix                 |
| -q         | Quality-score offset               |
| -Q         | Maximum quality score              |
| -n         | Maximum number of reads to process |
| -s         | Number of simulation iterations    |



### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Deployment

Add additional notes about how to deploy this on a live system

## Built With

* Python 3.7
* Numpy

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* **Darryn Zimire** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Stephens et al.,(2016) 



