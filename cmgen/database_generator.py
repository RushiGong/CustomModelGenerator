import yaml
from lxml import etree

def yaml_to_rng_input_strings(yaml_input_file):
    with open(yaml_input_file, 'r') as file:
        setting=yaml.safe_load(file)
    name=setting['database']["name"]
    Model_des=setting['database']["description"]
    Model_name=dict()
    Model_name[name]=Model_des
    Parameters=setting['database']["parameters"]
    if setting['database']["options"] is not None:
        Options=setting['database']["options"]
        return Model_name, Parameters, Options
    else: 
        return Model_name, Parameters
    
def generate_rng_schema(Model_name, Parameters, Options=None):
    model_name, model_des = list(Model_name.items())[0]
    # Define the RNG schema
    ns_structure = "http://relaxng.org/ns/structure/1.0"
    ns_annotations = "http://relaxng.org/ns/compatibility/annotations/1.0"
    nsmap = {None: ns_structure, 'a': ns_annotations}

    # Create grammar element
    grammar_schema = etree.Element("grammar", nsmap=nsmap)
    grammar_schema.set("datatypeLibrary", "http://www.w3.org/2001/XMLSchema-datatypes")    
    include_element = etree.SubElement(grammar_schema, "include", href="core.rng")

    #Model ConstituentArrary
    define_cons_element = etree.SubElement(grammar_schema, "define", name=model_name+"ConstituentArray")
    element_element = etree.SubElement(define_cons_element, 'element', name="ConstituentArray")
    one_or_more_element = etree.SubElement(element_element, 'oneOrMore')
    site_element = etree.SubElement(one_or_more_element, 'element', name="Site")
    
    id_attribute = etree.SubElement(site_element, 'optional')
    id_attr = etree.SubElement(id_attribute, 'attribute', name="id")
    id_data = etree.SubElement(id_attr, 'data', type="integer")
    
    ratio_attribute = etree.SubElement(site_element, 'optional')
    ratio_attr = etree.SubElement(ratio_attribute, 'attribute', name="ratio")
    ratio_data = etree.SubElement(ratio_attr, 'data', type="double")
    
    refid_attribute = etree.SubElement(site_element, 'optional')
    refid_attr = etree.SubElement(refid_attribute, 'attribute', name="refid")
    refid_data = etree.SubElement(refid_attr, 'data', type="integer")
    
    constituent_element = etree.SubElement(site_element, 'oneOrMore')
    constituent_sub_element = etree.SubElement(constituent_element, 'element', name="Constituent")
    constituent_attr = etree.SubElement(constituent_sub_element, 'attribute', name="refid")
    constituent_data = etree.SubElement(constituent_attr, 'data', type="string")
    
    #Model definition
    define_model_element = etree.SubElement(grammar_schema, "define", name=model_name+".model", combine="choice")
    model_element = etree.SubElement(define_model_element, 'element', name="Model")
    type_attribute = etree.SubElement(model_element, 'attribute', name="type")
    type_a_documentation = etree.SubElement(type_attribute, etree.QName(ns_annotations, "documentation"))
    type_a_documentation.text = model_des
    type_value = etree.SubElement(type_attribute, 'value')
    type_value.text = model_name

    interleave_element = etree.SubElement(model_element, 'interleave')
    model_constituent_array_ref = etree.SubElement(interleave_element, 'ref', name=model_name+"ConstituentArray")

    optional_element = etree.SubElement(interleave_element, 'optional')
    chemical_groups_element = etree.SubElement(optional_element, 'element', name="ChemicalGroups")
    chemical_groups_a_documentation = etree.SubElement(chemical_groups_element, etree.QName(ns_annotations, "documentation"))
    chemical_groups_a_documentation.text = "Mapping of species to integer chemical groups. Equal integers mean the species belong to the same group."

    one_or_more_element = etree.SubElement(chemical_groups_element, 'oneOrMore')
    constituent_element = etree.SubElement(one_or_more_element, 'element', name="Constituent")
    groupid_attr = etree.SubElement(constituent_element, 'attribute', name="groupid")
    groupid_data = etree.SubElement(groupid_attr, 'data', type="integer")
    refid_attr = etree.SubElement(constituent_element, 'attribute', name="refid")
    refid_data = etree.SubElement(refid_attr, 'data', type="string")
    
    #Parameter definition
    zero_or_more_element = etree.SubElement(define_model_element, 'zeroOrMore')
    parameter_element = etree.SubElement(zero_or_more_element, 'element', name="Parameter")
    type_attribute = etree.SubElement(parameter_element, 'attribute', name="type")
    type_choice = etree.SubElement(type_attribute, 'choice')

    value_list = list(Parameters.keys())
    a_documentation_list = list(Parameters.values())
    
    for value, a_documentation in zip(value_list, a_documentation_list):
        value_element = etree.SubElement(type_choice, 'value')
        value_element.text = value
        a_documentation_element = etree.SubElement(type_choice, etree.QName(ns_annotations, "documentation"))
        a_documentation_element.text = a_documentation

    interleave_element = etree.SubElement(parameter_element, 'interleave')
    optional_element = etree.SubElement(interleave_element, 'optional')
    text_element = etree.SubElement(optional_element, 'text')
    zero_or_more_element = etree.SubElement(interleave_element, 'zeroOrMore')
    ref_element = etree.SubElement(zero_or_more_element, 'ref', name="Interval")
    model_constituent_array_ref = etree.SubElement(interleave_element, 'ref', name=model_name+"ConstituentArray")

    
    #Optional addings
    if Options is not None:
        if len(Options) > 1 :
            choice_element = etree.SubElement(interleave_element, 'choice')
            for optadd in Options:
                optional_element = etree.SubElement(choice_element, 'optional')
                optadd_element = etree.SubElement(optional_element, 'element', name=optadd)
                comment=etree.Comment('Please finalize details of this optional adding.')
                optadd_element.append(comment)
        else:
            optadd = list(Options)[0]
            optional_element = etree.SubElement(interleave_element, 'optional')
            optadd_element = etree.SubElement(optional_element, 'element', name=optadd)
            comment=etree.Comment('Please finalize details of this optional adding.')
            optadd_element.append(comment)

    # Create the RNG schema tree
    comment=etree.Comment('Please modify database.rng and parser.py in pycalphad-xml')
    grammar_schema.append(comment)
    schema_tree = etree.ElementTree(grammar_schema)
    
    # Return the RNG schema tree
    return schema_tree

def save_rng_schema(schema_tree, filename):
    # Save the RNG schema to a file
    etree.indent(schema_tree, space="    ")
    with open(filename, "wb") as f:
        f.write(etree.tostring(schema_tree, encoding='utf-8', xml_declaration=True, pretty_print=True))
