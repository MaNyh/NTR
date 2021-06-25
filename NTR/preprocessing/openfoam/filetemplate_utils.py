import os



def readtxtfile(path_to_file):
    basepath = os.path.abspath(os.path.dirname(__file__))
    with open(os.path.join(basepath,path_to_file), "r") as fobj:
        content = fobj.readlines()
    return "".join(content)


def get_template_contents(templatepath, file_templates):
    template_contents = {}
    for key, vallist in file_templates.items():
        template_contents[key] = {}
        for val in vallist:
            template_contents[key][val] = readtxtfile(os.path.join(templatepath,key, val))

    return template_contents
