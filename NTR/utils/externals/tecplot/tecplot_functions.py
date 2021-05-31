import os


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


def openTecplotFile(path):

    values = []
    var = []
    zones = []

    zone_bool = -1
    with open(path, 'r') as f:
        for line in f:
            if line.startswith('VARIABLES'):
                line_string = line.replace('\n', '').split('=')[-1].split(',')
                for string in line_string:
                    var.append(string.replace('"', ''))

            if line.startswith('ZONE'):

                zone_bool = zone_bool + 1

                zones.append(line.replace('\n', '').split('=')[-1].replace('"', ''))
                list_01 = []
                for i in range(len(var)):
                    list_01.append([])

                values.append(list_01)

            if line.startswith('ZONE') or line.startswith('VARIABLES') or line.startswith('TITLE'):
                pass

            else:
                line_values = line.split()
                for i in range(len(line_values)):
                    values[zone_bool][i].append(float(line_values[i]))

    return values


def getTecplotFileVarNames(path):
    var = []

    with open(path, 'r') as f:
        for line in f:
            if line.startswith('VARIABLES'):
                line = line.replace(' ', '')
                line_string = line.replace('\n', '').split('=')[-1].split(',')

                for string in line_string:
                    var.append(string.replace('"', ''))

    return var


def getTecplotFileZoneNames(path):
    var = []

    with open(path, 'r') as f:
        for line in f:
            if line.startswith('ZONE') or line.startswith('Zone') or line.startswith('zone'):
                line = line.replace(' ', '')
                line_string = line.replace('\n', '').split('"')[1]

                var.append(line_string)

    return var
