# Fusion ART Retrieval Model by Patrick Tjahjadi
## Note
This README is for binary one-hot encoding, found in codes "professor-fusionART.ipynb" and "compl-professor-fusionART.ipynb".

## Aim
The aim of this project is to store and retrieve data from a Fusion Adaptive Resonance Theory
(ART) model, given a dataset. The accuracy of the model in retrieving data can then be measured
for varying levels of noise.

## Prerequisites and Setting Up
This model is written in Jupyter Notebook version 4.4.0 using Python version 3.7.1. The following software are required
to run the model:
* Python 3
* Jupyter Notebook 

The dataset used is "SCSE ProfProfile.csv".

Please ensure that the Fusion ART API "fusionART.py" and "ARTfunc.py" are present when running the code, since 
they provide the functions to build the Fusion ART model. These API are created by Dr. Budhitama Subagdja.

The code for data storage and retrieval using Fusion ART models can be found in "professor-fusionART.ipynb"
and "compl-professor-fusionART.ipynb". The former is to retrieve data without complement-coding, while the latter
uses complement-coding.

The following Python libraries are installed and used:
* pandas
* random
* numpy
* itertools
* sklearn
* csv

There is a batch file "research_query.bat" that is used to run the research keywords query for the Fusion ART 
model, which will output text files for varying levels of noise containing the results of the query.

## Running the code
Please run through the code in full to pre-process the dataset, build the Fusion ART model, store the data 
to the model and retrieve data for all noise levels.

The pre-processing involves converting the categorical data to binary data using one-hot encoding
to be stored to the model.

The following code is used to build the model:
```python
model = FusionART(schema = model_schema, beta = [1.0, 1.0, 1.0, 1.0], alpha = [0.1, 0.1, 0.1, 0.1], 
               gamma = [0.25, 0.25, 0.25, 0.25], rho = [0.2, 0.2, 0.5, 0.5])
model.F1Fields
```
The model parameters (alpha, gamma, rho) can be modified here for querying.

The data is stored to the model using the following code by Resonance Search:
```python
for i in range(0, len(name)):
    model.updateF1bySchema([{'name': data.columns[0], 'val': name_onehotlist[i]}, 
                            {'name': data.columns[1], 'val': group_onehotlist[i]}, 
                            {'name': data.columns[2], 'val': university_2d[i]}, 
                            {'name': data.columns[3], 'val': research_interest_2d[i]}])

    
    print("resonance search: ")
    J = model.resSearch()
    print("selected ", J)
    if model.uncommitted(J):
        print ('uncommitted')

    model.autoLearn(J)
```

There are two types of query retrievals to test the model. The first test involves retrieving the full data
of a professor, given their name as an input (query by name).
```python
def query_by_name(query_name)
```

The second query involves reading a set of research keyword interests, and retrieve the professor that
best matches the keywords despite noise (query by research keywords).
```python
def query_by_research_with_noise(query_research, output_file, noise_limit)
```
## Visuals
![Capture](https://user-images.githubusercontent.com/41354958/84756037-7ea35f00-afec-11ea-8cc8-f893d3c9c647.JPG)

Example of results of querying by name.

![Capture2](https://user-images.githubusercontent.com/41354958/84756126-98dd3d00-afec-11ea-82ed-dbf4500f37ba.JPG)

Example of results of querying by research keywords, for non-complement-coded data with 30% noise.

## Applications
We can retrieve the accuracy of the Fusion ART model in querying research keywords for varying levels of noise.
The model parameters affect the retrieval accuracy.

By modifying the pre-processing and querying process to account for different variables,
this model can be used for a variety of datasets.

## Special Thanks
* Prof. Tan Ah-Hwee, Singapore Management University
* Dr. Budhitama Subagdja, Nanyang Technological University
