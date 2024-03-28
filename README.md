# CustomModelGenerator

This **template generator** can help you create the templates to implement your custom model into [pycalphad](https://pycalphad.org/docs/latest/).

## How to start
1. Use the **[`Custom_Model_Database_Generator.ipynb`](./Custom_Model_Database_Generator.ipynb)** to create an  XML schema template for defining the thermodynamic database for your custom model.
2. Use the **[`Custom_Model_Template_Generator.ipynb`](./Custom_Model_Template_Generator.ipynb)** to create a PyCalphad-style model template for expressing your custom model in PyCalphad.

## Prepare configuration yaml file
See  **[`CustomModel.yaml`](./CustomModel.yaml)** for an example.

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
The ```model``` key is intended to define the thermodynamic model you want to implement
### name
The name of the thermodynamic model.

*type:* string

### energy_contributions
The ```energy_contributions``` keyword is intended to define what contributions to Gibbs energy are considered in the thermodynamic model

*type:* ```key:value``` pair<br>
```key``` could be the short name and ```value``` is the full name for the ```energy_function```

### basic_functions
Define some functions you would like to extract from existing models in pycalphad to use in the custom model.

*type:* list<br>
*default:* The minimum functions should be loaded from  **[`template_functions.json`](./template_functions.json)**

