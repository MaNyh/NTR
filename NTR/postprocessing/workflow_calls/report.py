import os
import itertools

from fpdf import FPDF

class PDF(FPDF):
    def __init__(self):
        super().__init__()
        self.WIDTH = 210
        self.HEIGHT = 297

    def header(self):
        # Custom logo and positioning
        # Create an `assets` folder and put any wide and short image inside
        # Name the image `logo.png`
        self.image('Logo_TFD_de.jpg', 10, 5, 33)
        self.set_font('Arial', 'B', 11)
        self.cell(self.WIDTH - 80)
        self.cell(60, 1, 'Sensivity Study', 0, 0, 'R')
        self.ln(20)

    def footer(self):
        # Page numbers in the footer
        self.set_y(-15)
        self.set_font('Arial', 'I', 8)
        self.set_text_color(128)
        self.cell(0, 10, 'Page ' + str(self.page_no()), 0, 0, 'C')

    def page_body(self, images):
        # Determine how many plots there are per page and set positions
        # and margins accordingly
        #if len(images) == 3:
        #    self.image(images[0], 15, 25, self.WIDTH - 60)
        #    self.image(images[1], 15, self.WIDTH / 2 + 5, self.WIDTH - 60)
        #    self.image(images[2], 15, self.WIDTH / 2 + 90, self.WIDTH - 60)
        if len(images) == 2:
            self.image(images[0], 15, 15, self.WIDTH - 60)
            self.image(images[1], 15, self.HEIGHT / 2 -15, self.WIDTH - 60)
        #else:
        #    self.image(images[0], 15, 25, self.WIDTH - 30)

    def print_page(self, images):
        # Generates the report
        self.add_page()
        self.page_body(images)


def construct(profile_pressure, entropydiff):

    # Construct data shown in document
    counter = 0
    pages_data = []
    temp = []
    # Get all plots
    a = [os.path.join(profile_pressure,i) for i in os.listdir(profile_pressure)]
    b = [os.path.join(entropydiff,i) for i in os.listdir(entropydiff)]

    files = list(itertools.chain.from_iterable(zip(a, b)))

    # Iterate over all created visualization
    for i in files:
        # We want 2 per page
        if counter == 2:
            pages_data.append(temp)
            temp = []
            counter = 0

        temp.append(i)
#        temp.append(b[i])
#        print(b[i])
        counter += 1

    return [*pages_data, temp]

def create_report():
    pdf = PDF()
    entropydiff = r"D:\Nachlauf_Sensitivitaeten\06_Contourplots"
    profile_pressure = r"D:\Nachlauf_Sensitivitaeten\05_Plots"
    plots_per_page = construct(profile_pressure, entropydiff)

    for elem in plots_per_page:
        pdf.print_page(elem)

    pdf.output('sensivity.pdf', 'F')

create_report()
