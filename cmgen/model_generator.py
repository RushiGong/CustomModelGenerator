import yaml
import os
import json

def custom_model(model_name, contributions, basic_functions, parameter_functions, energy_functions):
    # Create the class definition as a string
    class_definition = f"class {model_name}(Model):\n"
  
    # Add contributions to the class definition
    for contribution in contributions:
        class_definition += f"    {contribution}\n"

    # Add basic functions to the class definition
    for basic_function in basic_functions:
        class_definition += f"    {basic_function}\n"
   
    # Add parameter functions to the class definition
    for method in parameter_functions:
        class_definition += f"    {method}\n"
  
    for method in energy_functions:
        class_definition += f"    {method}\n"
    # Return the class definition
    return class_definition


def yaml_to_basic_functions_strings(setting):
    function_strings = []
    template_file=open('template_functions.json') ##To do: need to be more general.
    template_functions=json.load(template_file)
    for key in setting['model']['basic_functions']:
        if key == 'none':
            break
        if key == 'default':
            for fun_name in ['__new__', '__init__', '__eq__', '__ne__', '__hash__' 'moles', 'ast', 'variables', 'degree_of_ordering', 'quantities', 'endmember_reference_model', 'get_internal_constraints', '_array_validity', '_purity_test', '_interaction_test', '_site_ratio_normalization', 'redlich_kister_sum', 'build_phase']:
                for funs in template_functions['functions']:
                    if funs['name'] == fun_name: 
                        function_strings.append(funs['content']+'\n') 
        else:
            for funs in template_functions['functions']:
                if funs['name'] == key:
                    function_strings.append(funs['content']+'\n')
    return function_strings


def yaml_to_parameter_functions_strings(setting):
    function_strings = []
    functions=[]

    for param in setting['model']['parameters_functions']:
        if param['attributes'] is not None:
            def_string= f"def {param['parameter']}(self, dbe, {', '.join(param['attributes'])}):"
        else:
            def_string= f"def {param['parameter']}(self, dbe):"
        function_strings.append(def_string)
        if param['database_keyword'] is not None:
            keyword=param['database_keyword']
            search_string= f'\tparam_query=(\n\t\t\t(where("phase_name") == self.phase_name) & \\\n\t\t\t(where("parameter_type") == "{keyword}") & \\\n\t\t\t(where("constituent_array").test(self._array_validity))\n\t\t)\n\t\tparams = dbe._parameters.search(param_query)'
            function_strings.append(search_string)
        if param['comments'] is not None:
            comments_string=f"\t#{param['comments']}"
            function_strings.append(comments_string)
        function_strings.append(f"\treturn {param['parameter']}\n")
        functions.append(function_strings)
    
    return function_strings

def yaml_to_energy_functions_strings(setting):
    function_strings = []
    template_file=open('template_functions.json') ##To do: need to be more general.
    template_functions=json.load(template_file)
    for ene_f in setting['model']['energy_functions']:
        if ene_f['comments'] is not None:
            comments_string=ene_f['comments']
            function_strings.append(f'#{comments_string}\n')
        if ene_f['function'] == 'CEF-default':
            print('Need to include redlich_kister_sum function')
            for funs in template_functions['functions']:
                if funs['name'] == ene_f['energy']:
                    content_string=funs['content']
                    function_strings.append(f'{content_string}\n')
        else:
            function_strings.append(f"\n    def {ene_f['energy']}(self, dbe):\n")
            function_strings.append(f"\t{ene_f['function']}\n\t\treturn {ene_f['energy']}\n")
    return function_strings


def process_model_information(yamlfile):
    with open(yamlfile, 'r') as file:
        setting=yaml.safe_load(file)
    class_name=setting['model']['name']
    #contribution
    energy_contributions=setting['model']['energy_contributions']
    en_cons = [f'\t("{{key}}", "{{value}}")'.format(key=key, value=value) for key, value in energy_contributions.items()]
    energy_contributions_result_string = ["contributions = [\n{}\n\t]".format(",\n".join(en_cons))]
    #basic functions
    basic_function_strings=yaml_to_basic_functions_strings(setting)
    #parameter functions
    parameter_function_strings=yaml_to_parameter_functions_strings(setting)
    #energy functions
    energy_function_strings=yaml_to_energy_functions_strings(setting)
    return class_name, energy_contributions_result_string, basic_function_strings, parameter_function_strings, energy_function_strings

def model_generator(configuration_file, output_file, print_model=False):
    class_name, contributions, basic_functions, parameter_functions, energy_functions = process_model_information(configuration_file)
    dynamic_class_code = custom_model(class_name, contributions, basic_functions, parameter_functions, energy_functions)
    if print_model == True:
        print('Custom model template is:\n', dynamic_class_code)
    with open(output_file, "w") as py_file:
        f=open('./template_functions.json')
        temps=json.load(f)
        py_file.write(temps['imports']+'\n')
        py_file.write(dynamic_class_code)

