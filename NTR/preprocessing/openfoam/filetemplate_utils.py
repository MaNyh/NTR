import os

from NTR.utils.filehandling import readtxtfile


def get_template_contents(templatepath, file_templates):
    template_contents = {}
    for key, vallist in file_templates.items():
        template_contents[key] = {}
        for val in vallist:
            template_contents[key][val] = "".join(readtxtfile(os.path.join(templatepath, key, val)))

    return template_contents
