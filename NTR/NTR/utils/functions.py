import os
import yaml
import pickle



def yaml_dict_read(yml_file):

    args_from_yaml= {}

    with open(yml_file, "r") as Fobj:
        document = yaml.load_all(Fobj,Loader=yaml.FullLoader)
        for settings in document:
            for key, value in settings.items():
                args_from_yaml[key]=value

    return args_from_yaml

def write_igg_config(file, args):
    with open(file, "wb") as Fobj:
        pickle.dump(args, Fobj, protocol=0)


def writeTecplot1DFile(output_path, var_string, zone_string, values, title):
    # var_string: namen der variablen als liste ['U','p']
    # zone_string: namen der zonen als liste ['saugseute','druckseite']
    # values: erster index der liste steht fuer zone, dann folgen die listen der eigentlichen variablen
    # Beispiel: [[[10,11,10],[10000,11000,12000]],[[10,11,10],[10000,11000,12000]]]
    data = open(os.path.join(output_path), 'w')
    data.write('TITLE     ="' + title + '"\n')
    var = 'VARIABLES = '
    for i in range(len(var_string)):
        if i < len(var_string) - 1:
            var = var + '"' + var_string[i] + '", '
        else:
            var = var + '"' + var_string[i] + '"\n'

    data.write(var)

    for i in range(len(values)):
        data.write('ZONE T="' + zone_string[i] + '",I=' + str(len(values[i][0])) + '\n')

        # data.write('ZONE  T= "'+zone_string[i]+'"\n')
        for j in range(len(values[i][0])):
            line_string = ''
            for k in range(len(values[i])):
                if k < len(values[i]) - 1:
                    line_string = line_string + str(values[i][k][j]) + '\t'
                else:
                    line_string = line_string + str(values[i][k][j]) + '\n'
            data.write(line_string)

    data.close()


def writeTecplot2D3DFile(filename, X, Y, Z, vars):
    def pad(s, width):
        s2 = s
        while len(s2) < width:
            s2 = ' ' + s2
        if s2[0] != ' ':
            s2 = ' ' + s2
        if len(s2) > width:
            s2 = s2[:width]
        return s2

    def varline(vars, id, fw):
        s = ""
        for v in vars:
            s = s + pad(str(v[1][id]), fw)
        s = s + '\n'
        return s

    fw = 10  # field width

    f = open(filename, "wt")

    f.write('Variables="X","Y"')
    if len(Z) > 0:
        f.write(',"Z"')
    for v in vars:
        f.write(',"%s"' % v[0])
    f.write('\n\n')

    f.write('Zone I=' + pad(str(len(X)), 6) + ',J=' + pad(str(len(Y)), 6))
    if len(Z) > 0:
        f.write(',K=' + pad(str(len(Z)), 6))
    f.write(', F=POINT\n')

    if len(Z) > 0:
        id = 0
        for k in range(len(Z)):
            for j in range(len(Y)):
                for i in range(len(X)):
                    f.write(pad(str(X[i]), fw) + pad(str(Y[j]), fw) + pad(str(Z[k]), fw))
                    f.write(varline(vars, id, fw))
                    id = id + 1
    else:
        id = 0
        for j in range(len(Y)):
            for i in range(len(X)):
                f.write(pad(str(X[i]), fw) + pad(str(Y[j]), fw))
                f.write(varline(vars, id, fw))
                id = id + 1

    f.close()


def run_igg_meshfuncs():
    global args
    settings = yaml_dict_read("ressources/settings.yml")
    if settings["geom"]["create_bool"]:
        print("create_geometry")

        create(settings["geom"]["ptcloud_profile"],
               settings["geom"]["beta_meta_01"],
               settings["geom"]["beta_meta_02"],
               settings["geom"]["x_inlet"],
               settings["geom"]["x_outlet"],
               settings["geom"]["pitch"], )
    else:
        print("skipping geometry")
    if settings["mesh"]["create_bool"]:
        print("create_mesh")
        cwd = os.getcwd()
        os.chdir(settings["igg"]["install_directory"])
        igg_exe = settings["igg"]["executable"]
        script_path = os.path.join(cwd, "utils", "create_mesh.py")
        args_dict_path = os.path.join(cwd, settings["igg"]["argument_pickle_dict"])
        point_cloud_path = os.path.join(cwd, "ressources", "geom.dat")

        args = {}
        args["pointcloudfile"] = point_cloud_path
        write_igg_config(args_dict_path, args)

        os.system(igg_exe + " -batch -print -script " + script_path)
        os.remove(args_dict_path)
        # os.chdir(cwd)
    else:
        print("skipping meshing")
