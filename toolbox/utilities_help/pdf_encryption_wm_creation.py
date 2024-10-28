
from reportlab.pdfgen import canvas
from reportlab.lib.units import inch
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4


from pathlib import Path
from typing import Union, Literal, List
from PyPDF2 import PdfWriter, PdfReader

import argparse

# required input argument nome file input (file output default 'out.pdf')
# optional parameters (positional):  watermark (text of the watermark layer), output file name, print permission flag (T/F) edit permission flag (T/F), encryption password
# possibilità di creare un pdf/A


parser = argparse.ArgumentParser(description="required input argument: input file name. Optional input arguments (positional): watermark (text of the watermark layer), output file name, print permission flag (T/F) edit permission flag (T/F), encryption password.")
    
# Add an argument that accepts a variable number of parameters
parser.add_argument('parameters', nargs='+', help='Input parameters')
    
# Parse the arguments
args = parser.parse_args()

# debug code below
#for param in args.parameters:
#        print(param)
#print(len(args.parameters))

# default arguments
inputfile_text=args.parameters[0]
watermark_text="(C)DSconMATLAB"
outputfile_text="outputfile.pdf"
permissions=0b0000
password_text = 'DSconMATLAB24'

#if optional parameters are present

if len(args.parameters)==2:
    watermark_text=args.parameters[1]
   

if len(args.parameters)==3:
    watermark_text=args.parameters[1]
    outputfile_text=args.parameters[2]

if len(args.parameters)==4:
    watermark_text=args.parameters[1]
    outputfile_text=args.parameters[2]   
    print_flag=args.parameters[3]
    if print_flag=="T" or print_flag=="true":
        permissions=0b0100

if len(args.parameters)==5:
    watermark_text=args.parameters[1]
    outputfile_text=args.parameters[2]   
    print_flag=args.parameters[3]
    if print_flag=="T" or print_flag=="true":
        permissions=0b0100
    edit_flag=args.parameters[4]
    if edit_flag=="T" or edit_flag=="true":
        permissions=-1

if len(args.parameters)==6:
    watermark_text=args.parameters[1]
    outputfile_text=args.parameters[2]   
    print_flag=args.parameters[3]
    if print_flag=="T" or print_flag=="true":
        permissions=0b0100
    edit_flag=args.parameters[4]
    if edit_flag=="T" or edit_flag=="true":
        permissions=-1
    password_text=args.parameters[5]



def makeWatermark(text):
    pdf = canvas.Canvas("watermark_layer.pdf", pagesize=A4)
    pdf.translate(inch, inch)
    pdf.setFillColor(colors.grey, alpha=0.6)
    pdf.setFont("Helvetica", 85)
    pdf.rotate(45)
    pdf.drawCentredString(400, 100, text)
    pdf.save()



def watermark(
    content_pdf: Path,
    stamp_pdf: Path,
    pdf_result: Path,
    page_indices: Union[Literal["ALL"], List[int]] = "ALL",
):
    reader = PdfReader(content_pdf)
    if page_indices == "ALL":
        page_indices = list(range(0, len(reader.pages)))

    writer = PdfWriter()
    for index in page_indices:
        content_page = reader.pages[index]
        mediabox = content_page.mediabox

        # You need to load it again, as the last time it was overwritten
        reader_stamp = PdfReader(stamp_pdf)
        image_page = reader_stamp.pages[0]

        image_page.merge_page(content_page)
        image_page.mediabox = mediabox
        writer.add_page(image_page)

  
    # 0b0000 sets to 0 most bits of the permissions_flag de facto restricting any action on the PDF content.
    # note that the variable 'permissions_flag' is treated like a 16 bit binary flag table 
    # replicating Table 3.20 of the PDF 1.7 specification.
    writer.encrypt(user_password='', owner_pwd=password_text, permissions_flag=permissions)


    with open(pdf_result, "wb") as fp:
        writer.write(fp)


# begin code
makeWatermark(watermark_text)
watermark(inputfile_text, "watermark_layer.pdf", outputfile_text, "ALL")