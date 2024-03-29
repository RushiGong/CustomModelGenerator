# CustomModelGenerator

This **template generator** can help you create the templates to implement your custom thermodynamic model into [pycalphad](https://pycalphad.org/docs/latest/).

## How to start
1. Use the **[`Custom_Model_Database_Generator.ipynb`](./example/Custom_Model_Database_Generator.ipynb)** to create an  XML schema template for defining the thermodynamic database for your custom model.
2. Use the **[`Custom_Model_Template_Generator.ipynb`](./example/Custom_Model_Template_Generator.ipynb)** to create a PyCalphad-style model template for expressing your custom model in PyCalphad.

## Prepare configuration yaml file
See  **[`CustomModel.yaml`](./example/CustomModel.yaml)** for an example.

These are all key value pairs in the format:
```yaml
  key: value
```
They are nested for purely organizational purposes:
```yaml
  top_level_key:
    key: value
```
As long as keys are nested under the correct heading, they have no required order.
All of the possible keys are
```yaml
    model:
      name
      energy_contributions
      basic_functions
      parameters_functions
        - parameter
          attributes
          database_keyword
          comments
      energy_functions
        - energy
          function
          comments

    database:
      name
      description
      parameters
      options
```

## model
The ```model``` key is intended to define the custom model you want to implement.
### name
The name of the custom model.

*type:* string

### energy_contributions
The ```energy_contributions``` keyword is intended to define what contributions to Gibbs energy are considered in the custom model.

*type:* ```key:value``` pair<br>
```key``` could be the short name and ```value``` is the full name for the ```energy_function```.

### basic_functions
Define some functions you would like to extract from existing models in pycalphad to use in the custom model. See **[`template_functions.json`](./cmgen/template_functions.json)** for all available functions.

*type:* list<br>
*default:* The minimum functions should be loaded from  **[`template_functions.json`](./cmgen/template_functions.json)**.

### parameters_functions
The ```parameters_functions``` is intended to define the functions for new parameters in the custom model, you could provide information including parameter name, attributes, corresponding keyword defined in the database, and other comments for the parameter. 

#### parameter
Define the name of the new parameter.

*type:* string

#### attributes
Define the attributes of the parameter, such as this parameter is related to component i.

*type:* string

#### database_keyword
Match the parameter with the keyword defined in the database. 

*type:* string

#### comments
Additional information for the parameter.

*type:* string

### energy_functions
The ```energy_functions``` is intended to define the functions for energy contributions in the custom model, you could provide information including function name, function expression, and other comments.

#### energy
The ```energy``` keyword is corresponding with the ```value``` you defined in ```energy_contributions```, representing the full name for the ```energy_function```.

*type:* string

#### function
Enter the expression of energy function using this keyword. If ```CEF-default``` is used here, the same energy functions in the CEF model will be applied here. 

*type:* string

#### comments
Additional information for the energy function.

*type:* string

## database
The ```database``` key is intended to define the XML database schema you want to define for the custom model. 
### name
The name of the custom model. The same as the ```name``` value under the ```model``` key.

*type:* string

### description
Provide any comments or descriptions about your custom model, such as the full name of the model.

*type:* string

### parameters
Define the name of new parameters in the custom model. It is corresponding with the ```database_keyword``` value under ```model: parameters_functions``` key.

*type:* ```key:value``` pair<br>
```key``` could be the short name of the parameters and ```value``` is the full name or other comments for the parameters.

### options
Define some optional schema for parameter structures. Options include ```Exponent``` and ```Order```.

*type:* string
