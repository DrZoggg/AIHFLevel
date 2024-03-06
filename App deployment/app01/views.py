import subprocess
import time
from django.shortcuts import render,redirect,HttpResponse
from app01 import models
from django import forms
from django.core.validators import MinValueValidator, MaxValueValidator,RegexValidator




class MyForm(forms.Form):
    age = forms.IntegerField(label='Age', validators=[MinValueValidator(18), MaxValueValidator(95)])
    arrhythmia_chioces = (
        (1,"yes"),
        (0, "no"),
    )
    arrhythmia= forms.ChoiceField(label="Co-morbidity(Arrhythmia)",choices=arrhythmia_chioces)

    CAD_chioces = (
        (1,"yes"),
        (0, "no"),
    )


    CAD = forms.ChoiceField(label="Co-morbidity(CHD)",choices=CAD_chioces)

    CKD_stage_chioces= (
        (1,"I Stage"),
        (2, "II Stage"),
        (3, "III Stage"),
        (4, "IV Stage"),
        (5, "V Stage"),
    )


    CKD_stage = forms.ChoiceField(label="CKD stage",choices=CKD_stage_chioces,initial=1)

    Cr = forms.FloatField(label="Cr",validators=[MinValueValidator(25.0), MaxValueValidator(2200.0)])
    EF = forms.FloatField(label="LVEF",validators=[MinValueValidator(10.0), MaxValueValidator(95.0)])
    GFR = forms.FloatField(label="eGFR",validators=[MinValueValidator(5.000), MaxValueValidator(130.000)])
    LY_Per = forms.FloatField(label="Lymphocyte",validators=[MinValueValidator(1.0), MaxValueValidator(100.0)])
    MCHC = forms.FloatField(label="MCHC",validators=[MinValueValidator(250.000), MaxValueValidator(400.000)])
    SV = forms.FloatField(label="SV",validators=[MinValueValidator(10.0), MaxValueValidator(250.0)])
    TN = forms.FloatField(label="CTnI",validators=[MinValueValidator(0.00), MaxValueValidator(50.00)])
    TBIL = forms.FloatField(label="TBIL",validators=[MinValueValidator(0.50), MaxValueValidator(100.00)])
    Timepoint = forms.IntegerField(label="Indicated time point",validators=[MinValueValidator(1), MaxValueValidator(939)])


    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        for name, field in self.fields.items():
            if field.label == "Age":
                field.widget.attrs = {"class": "form-control", "placeholder": "Enter 18 - 95"}
            elif field.label == "Cr":
                field.widget.attrs = {"class": "form-control", "placeholder": "Enter 25.0 - 2200.0"}
            elif field.label == "LVEF":
                field.widget.attrs = {"class": "form-control", "placeholder": "Enter 10.0 - 95.0"}
            elif field.label == "eGFR":
                field.widget.attrs = {"class": "form-control", "placeholder": "Enter 5.000 - 130.000"}
            elif field.label == "Lymphocyte":
                field.widget.attrs = {"class": "form-control", "placeholder": "Enter 1.0 - 100.0"}
            elif field.label == "MCHC":
                field.widget.attrs = {"class": "form-control", "placeholder": "Enter 250.0 - 400.0"}
            elif field.label == "SV":
                field.widget.attrs = {"class": "form-control", "placeholder": "Enter 10.0 - 250.0"}
            elif field.label == "CTnI":
                field.widget.attrs = {"class": "form-control", "placeholder": "Enter 0.00 - 50.00"}
            elif field.label == "TBIL":
                field.widget.attrs = {"class": "form-control", "placeholder": "Enter 0.50 - 100.0"}
            elif field.label == "Indicated time point":
                field.widget.attrs = {"class": "form-control", "placeholder": "Enter 1 - 939"}
            else:
                field.widget.attrs = {"class": "form-control", "placeholder": "Enter"}





def runcmd(cmd):
    result = []
    process = subprocess.Popen(cmd,
                               stdin=subprocess.PIPE,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               # env=my_env,
                               shell=True,
                               encoding='UTF8')
    for line in process.stdout:
        result.append(line)

    output, errcode = process.communicate()
    if errcode == 0:
        print('error code is 0')
    else:
        print('success---!', errcode)
        print(output)
    return result


def submit(request):
    if request.method == "GET":
        form = MyForm()
        return render(request, 'Submit.html', {"form": form})
    form = MyForm(data=request.POST)
    if form.is_valid():
        print("进var：",time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time())))
        print("555555",form.cleaned_data["age"],form.cleaned_data["arrhythmia"],form.cleaned_data["CAD"],form.cleaned_data["CKD_stage"],form.cleaned_data["Cr"],form.cleaned_data["EF"],form.cleaned_data["GFR"],form.cleaned_data["LY_Per"],
              form.cleaned_data["MCHC"],form.cleaned_data["SV"],form.cleaned_data["TN"],form.cleaned_data["TBIL"],form.cleaned_data["Timepoint"])
        dir_result="/static/R/result_img/"
        filename=str(form.cleaned_data["age"])+str(form.cleaned_data["arrhythmia"])+str(form.cleaned_data["CAD"])+str(form.cleaned_data["CKD_stage"])\
               +str(form.cleaned_data["Cr"])+str(form.cleaned_data["EF"])+str(form.cleaned_data["GFR"])+str(form.cleaned_data["LY_Per"])+str(form.cleaned_data["MCHC"])\
               +str(form.cleaned_data["SV"])+str(form.cleaned_data["TN"])+str(form.cleaned_data["TBIL"])+str(form.cleaned_data["Timepoint"])
        filename1= dir_result+filename+"_imp_radar.png"
        filename2= dir_result+filename+"_imp_bar.png"
        filename3= dir_result+filename+"_KMplot.png"
        filename4= dir_result+filename+"_Survplot.png"

        char = "Rscript /home/web/Rnamespace/Web_Rcode1127.R -A "+str(form.cleaned_data["age"])+" -H "+str(form.cleaned_data["arrhythmia"])+" -C "+str(form.cleaned_data["CAD"])+" -K "+str(form.cleaned_data["CKD_stage"])\
               +" -R "+str(form.cleaned_data["Cr"])+" -E "+str(form.cleaned_data["EF"])+" -G "+str(form.cleaned_data["GFR"])+" -L "+str(form.cleaned_data["LY_Per"])+" -M "+str(form.cleaned_data["MCHC"])+" -S "\
               +str(form.cleaned_data["SV"])+" -T "+str(form.cleaned_data["TN"])+" -B "+str(form.cleaned_data["TBIL"])+" -I "+str(form.cleaned_data["Timepoint"])+" -F "+filename
        res = runcmd(char)
        res1=res[-5:]
        risk = res1[0]
        risk = risk[5:-2]
        half = res1[1]
        half = half[4:]
        year1= res1[2]
        year1 = year1[4:]
        year2 = res1[3]
        year2 = year2[4:]
        day900 = res1[4]
        day900 = day900[4:]
        return render(request, 'result.html', {"form": form,"filename1":filename1,"filename2":filename2,"filename3":filename3,"filename4":filename4
                                               ,"risk":risk,"half":half,"year1":year1,"year2":year2,"day900":day900})
    else:
        return render(request, 'Submit.html', {"form": form})



def result(request):
    if request.method == "GET":
     return render(request, 'result.html')



