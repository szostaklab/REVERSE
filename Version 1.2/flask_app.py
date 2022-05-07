
# coding:utf-8
from flask import Flask, request, make_response, render_template, redirect, url_for, session, flash, send_file, send_from_directory, jsonify, Response
from processing import do_calculation, start, end, cutoff_conserved, process_data
import io
import collections
from helpers import error_page
import numpy as np
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
from matplotlib.pyplot import figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
plt.rcParams["figure.figsize"] = [7.0, 3.50]
plt.rcParams["figure.autolayout"] = True
import matplotlib.pyplot as plt
import plotly.express as px
import seaborn as sns
import difflib
import os
import uuid
import shutil
import zipfile
from Bio.Seq import Seq
import matplotlib.gridspec as gridspec
import random
random.seed(10)
import re
import mechanize
from selenium import webdriver
from time import sleep
import http.cookiejar
import requests
from selenium.webdriver.chrome.options import Options
from waitress import serve
from selenium.webdriver.support.select import Select
#!/usr/bin/python2.7
import numpy as np
import scipy
import matplotlib.pyplot as plt
from sklearn.metrics import pairwise_distances #jaccard diss.
from sklearn import manifold  # multidimensional scaling
from mpl_toolkits.mplot3d import axes3d
from matplotlib import style
import glob
from multiprocessing import Lock
from multiprocessing.managers import AcquirerProxy, BaseManager, DictProxy
from collections import OrderedDict
import logging
from flask import Response
from flask import stream_with_context
import requests
from flask import request
from flask import render_template, Blueprint, request, make_response
from werkzeug.utils import secure_filename
from flask import Flask,render_template,request,redirect,url_for
from werkzeug.utils import secure_filename
import os
from time import sleep
from flask import copy_current_request_context
import threading
import datetime
import time
import statistics
import logomaker as lm
import pandas as pd
import seaborn
import matplotlib.pyplot as plt
plt.style.use('seaborn-ticks')
from matplotlib import transforms
import matplotlib.patheffects


app = Flask(__name__)
@app.route('/index', methods=['POST','GET'])
def index():
    @copy_current_request_context
    def save_file(closeAfterWrite):
        print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S') + " i am doing")
        f = request.files['file']
        basepath = os.path.dirname(__file__)
        upload_path = os.path.join(basepath, '',secure_filename(f.filename))
        f.save(upload_path)
        closeAfterWrite()
        print(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S') + " write done")
    def passExit():
        pass
    if request.method == 'POST':
        f= request.files['file']
        normalExit = f.stream.close
        f.stream.close = passExit
        t = threading.Thread(target=save_file,args=(normalExit,))
        t.start()
        return redirect(url_for('index'))
    return render_template('index.html')

def get_clusters(clustering):
    num_clustering = clustering
    seqs_use = shared_dict["all_high_qual_seqs"]
    shared_dict["all_labels"]=[ [] for _ in range(len(seqs_use)) ]
    all_labels = []
    for round_num in range(len(seqs_use)):
                high_quality_sequences = list(set(seqs_use[round_num]))
                matrix = np.asarray([np.fromstring(str(s), dtype=np.uint8) for s in high_quality_sequences]);
                error = 10e10
                lowest_k = 10e10
                for k in range(1, 1+num_clustering):
                  print('Tried k='+str(k))
                  kmeans = KMeans(init="random", n_clusters=k,n_init=10,max_iter=300,random_state=42)
                  kmeans.fit(matrix)
                  kmeans.labels_[:]
                  if kmeans.inertia_<error:
                    error = kmeans.inertia_
                    lowest_k = k
                kmeans = KMeans(init="random", n_clusters=lowest_k,n_init=10,max_iter=300,random_state=42)
                kmeans.fit(matrix)
                labels = kmeans.labels_[:]
                all_labels.append(list(labels))

    shared_dict["all_labels"] = all_labels

    return all_labels

def get_shared_state(host, port, key):
    shared_dict = {}
    shared_lock = Lock()
    manager = BaseManager((host, port), key)
    manager.register("get_dict", lambda: shared_dict, DictProxy)
    manager.register("get_lock", lambda: shared_lock, AcquirerProxy)
    try:
        manager.get_server()
        manager.start()
    except OSError:  # Address already in use
        manager.connect()
    return manager.get_dict(), manager.get_lock()

HOST = "127.0.0.1"
PORT = 35791
KEY = b"secret"

shared_dict, shared_lock = get_shared_state(HOST, PORT, KEY)


app = Flask(__name__)
app.config["DEBUG"] = True
app.config["SECRET_KEY"] = "lkmaslkdsldsamdlsdmasldsmkdd"

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
directory=BASE_DIR+'/files'


def get_uploaded_files(file):
    global uploaded_files
    os.chdir("/home/zweiss124/mysite/files")
    #for file in shared_dict["filenames"]:
        #uploaded_files.append(list(open("/home/zweiss124/mysite/files"+'/'+str(file)).readlines()[:100] ))
    uploaded_files = (list(open("/home/zweiss124/mysite/files"+'/'+str(file)).readlines()))
        #uploaded_files.append(list(open("/home/zweiss124/mysite/files"+'/'+str(file)).readlines()))
    return uploaded_files


def delete_files():
    filelist = glob.glob(os.path.join("/home/zweiss124/mysite/files", "*fastq"))
    for f in filelist:
        os.remove(f)
    filelist = glob.glob(os.path.join("/home/zweiss124/mysite/files", "*zip"))
    for f in filelist:
        os.remove(f)
    dir = "/home/zweiss124/mysite/files"
    for f in os.listdir(dir):
        if str(f[0]) != '_':
            os.remove(os.path.join(dir, f))


    os.chdir(directory)
    files=glob.glob('*.fastq')
    for filename in files:
        os.unlink(filename)

def get_high_qual_seqs(trim_type, filter,uploaded_files, start, end, rc, enter_qual_cutoff, rounds):

    all_seqs = []
    all_high_qual_seqs=[]


    for file in range(shared_dict["num_rounds"]):

        quality = uploaded_files[file][3::4]
        seqs = uploaded_files[file][1::4]

        ##UNUSED
        #Remove duplicates:
        unique_seqs = []
        unique_qualities = []
        indicies_of_unique_seqs = list(np.unique(seqs, return_index=True)[1])
        for index in sorted(np.unique(seqs, return_index = True)[1])[:-1]:
          unique_qualities.append(quality[index])
          unique_seqs.append(seqs[index])
        ##UNUSED

        all_seqs.append(seqs)

        origonal_seqs = all_seqs

        #NO FILTER BY QUALITY
        if filter == 'unfiltered':
            all_high_qual_seqs = all_seqs

        #FILTER BY QUALITY
        else:
            high_quality = []
            acceptable_quals_1_percent_error = [5,6,7,8,9,':',';', '<', '=', '>', '?', '@', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K']
            for seq in range(len(quality)):
              overlap=0
              for char in list(set(quality[seq]).intersection(acceptable_quals_1_percent_error)):
                    overlap+=quality[seq].count(char)

              if overlap/len(quality[1]) > shared_dict["quality_cutoff"]/100:
                    high_quality.append(seqs[seq].strip("\n"))
            all_high_qual_seqs.append(high_quality)



    if trim_type == 'pos':
        #TRIM
        all_trimmed_seqs = []
        for round in range(len(all_high_qual_seqs)):
            all_trimmed_seqs.append([w[start:end] for w in all_high_qual_seqs[round]])
        all_high_qual_seqs = all_trimmed_seqs

        orig_trimmed_seqs = []
        for round in range(len(origonal_seqs)):
            orig_trimmed_seqs.append([w[start:end] for w in origonal_seqs[round]])
        origonal_seqs = orig_trimmed_seqs

    if trim_type == 'seq':
        #TRIM
        all_trimmed_seqs = []
        for round in range(len(all_high_qual_seqs)):
            trimmed_seqs = []
            for w in all_high_qual_seqs[round]:
                if start in w and end in w:
                    trimmed_seqs.append(w[w.index(start)+len(start):w.index(end)])
            all_lens = [len(i) for i in trimmed_seqs]
            correct_len=[]
            for seq in trimmed_seqs:
                if len(seq)==statistics.mode(all_lens):
                    correct_len.append(seq)
            all_trimmed_seqs.append(correct_len)
            #all_trimmed_seqs.append(trimmed_seqs)
        all_high_qual_seqs = all_trimmed_seqs

        orig_trimmed_seqs = []
        for round in range(len(origonal_seqs)):
            trimmed_seqs = []
            for w in origonal_seqs[round]:
                if start in w and end in w:
                    trimmed_seqs.append(w[w.index(start)+len(start):w.index(end)])
            #all_lens = [len(i) for i in trimmed_seqs]
            #correct_len=[]
            #for seq in trimmed_seqs:
            #    if len(seq)==statistics.mode(all_lens):
            #        correct_len.append(seq)
            #orig_trimmed_seqs.append(correct_len)
            orig_trimmed_seqs.append(trimmed_seqs)
        origonal_seqs = orig_trimmed_seqs

   #RC
    if rc == 'convert':
        all_rc_seqs = []
        for round in range(len(all_high_qual_seqs)):
            all_rc_seqs.append([str(Seq(w).reverse_complement()) for w in all_high_qual_seqs[round]])
        all_high_qual_seqs = all_rc_seqs

        orig_rc_seqs = []
        for round in range(len(origonal_seqs)):
            orig_rc_seqs.append([str(Seq(w).reverse_complement()) for w in origonal_seqs[round]])
        origonal_seqs = orig_rc_seqs

    return all_high_qual_seqs[:10000], origonal_seqs[:10000]

def get_high_qual_seqs_NOTUNIQUE(trim_type, filter,uploaded_files, start, end, rc, enter_qual_cutoff, rounds):
    return get_high_qual_seqs(trim_type, filter,uploaded_files, start, end, rc, enter_qual_cutoff, rounds)


@app.route('/choose_rounds', methods=["GET", "POST"])
def choose_rounds():
    if request.method == 'POST':

        choose_rounds_start = time.time()

        delete_files()

        shared_dict["num_rounds"] = 1
        shared_dict["option"] = 'option1'
        shared_dict["num_end"] = 1000000
        shared_dict["num_start"] = 0
        shared_dict["seq_start"] = ''
        shared_dict["seq_end"] = ''
        shared_dict["rc"] = 0
        shared_dict["quality_cutoff"] = 0
        shared_dict["uploaded_files"]=[]
        shared_dict["all_high_qual_seqs"]=0
        shared_dict["origonal_seqs"]=0
        shared_dict["all_tracking_nums"]=[]
        shared_dict["all_tracking_counts"] = []
        shared_dict["filenames"] = ''
        shared_dict["num_to_upload"] = 100
        shared_dict["all_labels"] = []
        shared_dict["top_clusters"] = []
        shared_dict["all_diffs"] = []
        shared_dict["heatmap_data"] = []
        shared_dict["size_to_upload"] = ''
        shared_dict["round_num"] = 0
        shared_dict["counter"] = 0
        shared_dict["times"] = ''
        shared_dict["all_regions_all_counts_rounds_seqs"] = []
        shared_dict["qual_yes"] = ''
        shared_dict["top"] = []
        shared_dict["num_clusters"] = 0
        shared_dict["trim_type"] = 'pos'
        shared_dict["seq_start"] = ''
        shared_dict["seq_end"] = ''
        shared_dict["most_abundant_seqs"] = []
        shared_dict["dfObj"] = []
        shared_dict["all_scores"] = []

        if not request.form.get("numrds"):
            return error_page("you need to enter a number of rounds")

        if request.form.get("numrds").isdigit()==False:
            return error_page("you need to enter a valid number of rounds")

        if not request.form.get("qual_yes"):
            return error_page("you need to enter a quality track")

        shared_dict["qual_yes"] = request.form.get("qual_yes")

        if shared_dict["qual_yes"] == 'filtered':
            if not request.form.get("quality"):
                return error_page("you need to enter a quality cutoff")
            if request.form.get("quality").isdigit()==False:
                return error_page("you need to enter a valid quality cutoff")
            shared_dict["quality_cutoff"] = int(request.values.get("quality"))
            if type(shared_dict["quality_cutoff"])!=int:
                return error_page("you need to enter a valid quality cutoff")
            if shared_dict["quality_cutoff"]>100:
                return error_page("quality cutoff must be between 0 and 100")

        if not request.form.get("size_to_analyze"):
            return error_page("you need to choose a size of sequences to analyze")

        if not request.form.get("num_to_analyze"):
            return error_page("you need to choose a number of sequences to analyze")

        shared_dict["size_to_upload"] = request.form.get("size_to_analyze")

        shared_dict["num_to_upload"] = 4*int(request.form.get("num_to_analyze"))

        shared_dict["num_rounds"] = int(request.values.get("numrds"))

        option = 'option1'

        #quality_cutoff = 90



        if type(shared_dict["num_rounds"])!=int:
            return error_page("you need to enter a valid number of rounds")


        if shared_dict["num_rounds"]>10:
            return error_page("round number must be between 0 and 10")

        global counter
        counter = 0

        if shared_dict["size_to_upload"] == 'small':
            choose_rounds_end = time.time()
            shared_dict["times"]+='\n\Choose Rounds('
            shared_dict["times"]+=str(round(choose_rounds_end-choose_rounds_start, 3))
            shared_dict["times"]+=')'
            times = shared_dict["times"]
            return render_template("upload_files.html", times=times, size = shared_dict["size_to_upload"], num_to_analyze = shared_dict["num_to_upload"], num_rounds = shared_dict["num_rounds"], counter=counter, quality_cutoff=shared_dict["quality_cutoff"], option = shared_dict["option"])

        else:
            if shared_dict["size_to_upload"] == 'large':
                choose_rounds_end = time.time()
                shared_dict["times"]+=str('\n\\Choose Rounds(')
                shared_dict["times"]+=str(round(choose_rounds_end-choose_rounds_start, 3))
                shared_dict["times"]+=')'
                times = shared_dict["times"]
                return render_template("upload_files_LARGE.html", times = times, round_num=0, num_parts = 2, size = shared_dict["size_to_upload"], num_to_analyze = shared_dict["num_to_upload"], num_rounds = shared_dict["num_rounds"], counter=counter, quality_cutoff=shared_dict["quality_cutoff"], option = shared_dict["option"])
            else:
                return error_page("error in size of file")

    if request.method == "GET":
        return render_template("choose_rounds.html")

def unzip_file(zip_src, dst_dir):
    """
    Unzip the zip file
         :param zip_src: full path of zip file
         :param dst_dir: the destination folder to extract to
    :return:
    """
    r = zipfile.is_zipfile(zip_src)
    if r:
        fz = zipfile.ZipFile(zip_src, "r")
        for file in fz.namelist():
            fz.extract(file, dst_dir)
    else:
                 return "Please upload zip file"

@app.route("/upload", methods=["GET", "POST"])
def upload():
    if request.method == "GET":
        times = shared_dict["times"]
        return render_template("upload_files.html",times=times, num_rounds = shared_dict["num_rounds"],  counter=counter, quality_cutoff=shared_dict["quality_cutoff"], rounds=range(shared_dict["num_rounds"]), option = shared_dict["option"])

    upload_start = time.time()

    upload_start = time.time()
    global all_obj
    global all_target_paths
    global fname
    global uploaded_files

    if shared_dict["size_to_upload"] == 'large':

            all_target_paths = []
            all_obj = []


            os.chdir("/home/zweiss124/mysite/files")

            #obj = request.files.get("file_"+str(shared_dict["round_num"]))
            obj = request.files.get("file_0")

            ret_list = obj.filename.rsplit(".", maxsplit=1)
            file_path = os.path.join(BASE_DIR, "files", obj.filename) # The path where the uploaded file is saved
            fname = obj.filename

            obj.save(file_path)
            shared_dict["counter"]+=1

            if shared_dict["counter"]%2!=0:
                round_num=shared_dict["round_num"]
                return render_template("upload_files_LARGE.html", round_num=round_num, all_fnames = shared_dict["filenames"], num_parts = 2, size = shared_dict["size_to_upload"], num_to_analyze = shared_dict["num_to_upload"], num_rounds = shared_dict["num_rounds"], counter=shared_dict["counter"], quality_cutoff=shared_dict["quality_cutoff"], option = shared_dict["option"])

            if shared_dict["counter"]%2==0:
                os.system('cat round'+str(shared_dict["round_num"]+1)+'aa '+'round'+str(shared_dict["round_num"]+1)+"ab "+'>> combined_data.zip')
                file_path = os.path.join(BASE_DIR, "files", "combined_data.zip")
                target_path = os.path.join(BASE_DIR, "files") # The path where the unzipped files are saved
                ret = unzip_file(file_path, target_path)
                os.remove(file_path) # delete file

            for file in glob.glob("*.fastq"):
                shared_dict["filenames"]+=(str(file).strip(".zip"))+","

            files_on=0
            os.chdir("/home/zweiss124/mysite/files")
            for file in glob.glob("*.fastq"):
                files_on+=1
            global uf
            global test
            test = 0
            uf = shared_dict["uploaded_files"]


            if files_on == shared_dict["num_rounds"]:
                test = 1
                uploaded_files = []
                for file in range(shared_dict["num_rounds"]):
                    uploaded_files+=[get_uploaded_files(shared_dict["filenames"].split(',')[file])[:shared_dict["num_to_upload"] ]]



                shared_dict["uploaded_files"] = uploaded_files

                if len(shared_dict["uploaded_files"]) not in range(10):
                    return error_page("error uploading files - number of uploaded files not between 1-10")

            upload_end = time.time()

            shared_dict["times"]+=str('\\Upload(')
            shared_dict["times"]+=str(round(upload_end-upload_start, 3))
            shared_dict["times"]+=')'

            if shared_dict["num_rounds"]>1:
                if files_on != shared_dict["num_rounds"]:
                    shared_dict["counter"]=0
                    os.chdir("/home/zweiss124/mysite/files")
                    for file in glob.glob("*"):
                        shared_dict["counter"]+=1
                    upload_end = time.time()
                    shared_dict["times"]+=('\\Upload('+str(round(upload_end-upload_start, 3)))
                    times = shared_dict["times"]
                    shared_dict["times"]+=')'
                    shared_dict["round_num"] +=1
                    round_num = shared_dict["round_num"]
                    return render_template("upload_files_LARGE.html",round_num=round_num,times=times, test=test, uf = shared_dict["uploaded_files"], all_fnames = shared_dict["filenames"], filenames = fname, num_rounds=shared_dict["num_rounds"], counter=shared_dict["counter"], quality_cutoff=shared_dict["quality_cutoff"], option = shared_dict["option"])
                else:
                    return error_page(str(shared_dict["uploaded_files"]))
                    shared_dict["all_high_qual_seqs"], shared_dict["origonal_seqs"] = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], 0, 10000, 0, shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))
                    if len(shared_dict["all_high_qual_seqs"][0])==0:
                        return error_page("please choose a lower quality threshold")
                    upload_end = time.time()
                    shared_dict["times"]+=('\\Upload('+str(round(upload_end-upload_start, 3)))
                    shared_dict["times"]+=')'
                    times = shared_dict["times"]
                    return render_template("preview.html", times=times, num_rounds = shared_dict["num_rounds"], rounds=range(shared_dict["num_rounds"]), uploaded_files=shared_dict["uploaded_files"], all_high_qual_seqs=shared_dict["all_high_qual_seqs"], origonal_seqs=shared_dict["origonal_seqs"], quality_cutoff=shared_dict["quality_cutoff"] )

            else:

                shared_dict["all_high_qual_seqs"], shared_dict["origonal_seqs"] = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], 0, 10000, 0, shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))
                if len(shared_dict["all_high_qual_seqs"][0])==0:
                        return error_page("please choose a lower quality threshold")
                upload_end = time.time()
                shared_dict["times"]+=('\\Upload('+str(round(upload_end-upload_start, 3)))
                shared_dict["times"]+=')'
                times = shared_dict["times"]
                return render_template("preview.html", times=times, num_rounds = shared_dict["num_rounds"], rounds=range(shared_dict["num_rounds"]), uploaded_files=shared_dict["uploaded_files"], all_high_qual_seqs=shared_dict["all_high_qual_seqs"], origonal_seqs=shared_dict["origonal_seqs"], quality_cutoff=shared_dict["quality_cutoff"] )



    if shared_dict["size_to_upload"] != 'small' and shared_dict["size_to_upload"] !='large':
        return error_page("unknown file size")

    if shared_dict["size_to_upload"] == 'small':

            all_target_paths = []
            all_obj = []

            shared_dict["round_num"] = 0

            obj = request.files.get("file_"+str(shared_dict["round_num"]))
            ret_list = obj.filename.rsplit(".", maxsplit=1)
            if len(ret_list) != 2:
                            return error_page("please upload zip file")
            if ret_list[1] != "zip":
                            return error_page("please upload zip file")

            file_path = os.path.join(BASE_DIR, "files", obj.filename) # The path where the uploaded file is saved

            fname = obj.filename

            shared_dict["filenames"]+=(str(fname).strip(".zip"))+","

            obj.save(file_path)
            target_path = os.path.join(BASE_DIR, "files") # The path where the unzipped files are saved
            ret = unzip_file(file_path, target_path)
            os.remove(file_path) # delete file


            files_on=0
            os.chdir("/home/zweiss124/mysite/files")
            for file in glob.glob("*.fastq"):
                files_on+=1


            test = 0
            uf = shared_dict["uploaded_files"]
            if files_on == shared_dict["num_rounds"]:
                test = 1

                uploaded_files = []

                for file in range(shared_dict["num_rounds"]):
                    uploaded_files+=[get_uploaded_files(shared_dict["filenames"].split(',')[file])[:shared_dict["num_to_upload"] ]]

                shared_dict["uploaded_files"] = uploaded_files

                if len(shared_dict["uploaded_files"]) not in range(10):
                    return error_page("error uploading files - number of uploaded files not between 1-10")

            if shared_dict["num_rounds"]>1:
                shared_dict["counter"]=0
                os.chdir("/home/zweiss124/mysite/files")
                for file in glob.glob("*.fastq"):
                    shared_dict["counter"]+=1
                upload_end = time.time()
                shared_dict["times"]+=('\\Upload('+str(round(upload_end-upload_start, 3)))
                shared_dict["times"]+=')'
                times = shared_dict["times"]
                return render_template("upload_files.html",times=times, test=test, uf = shared_dict["uploaded_files"], all_fnames = shared_dict["filenames"], filenames = fname, num_rounds=shared_dict["num_rounds"], counter=shared_dict["counter"], quality_cutoff=shared_dict["quality_cutoff"], option = shared_dict["option"])

            else:
                shared_dict["all_high_qual_seqs"], shared_dict["origonal_seqs"] = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], 0, 10000, 0, shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))
                if len(shared_dict["all_high_qual_seqs"][0])==0:
                        return error_page("please choose a lower quality threshold")
                upload_end = time.time()
                shared_dict["times"]+=('\\Upload('+str(round(upload_end-upload_start, 3)))
                shared_dict["times"]+=')'
                times = shared_dict["times"]
                return render_template("preview.html", times=times, num_rounds = shared_dict["num_rounds"], rounds=range(shared_dict["num_rounds"]), uploaded_files=shared_dict["uploaded_files"], all_high_qual_seqs=shared_dict["all_high_qual_seqs"], origonal_seqs=shared_dict["origonal_seqs"], quality_cutoff=shared_dict["quality_cutoff"] )

round_counter = 0

@app.route('/finalize', methods=["GET", "POST"])
def finalize():
    if request.method == 'GET':
        times = shared_dict["times"]
        return render_template("upload_files.html", times=times, num_rounds=shared_dict["num_rounds"], counter=shared_dict["counter"], quality_cutoff=shared_dict["quality_cutoff"], rounds=range(shared_dict["num_rounds"]), option = shared_dict["option"])

    if request.method == 'POST':
        qual_start = time.time()
        shared_dict["all_high_qual_seqs"], shared_dict["origonal_seqs"] = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], 0, 10000, 0, shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))
        qual_end = time.time()

        if len(shared_dict["all_high_qual_seqs"][0])==0:
                        return error_page("please choose a lower quality threshold")
        rounds = range(shared_dict["num_rounds"])

        global len_uploaded_files
        len_uploaded_files = len(shared_dict["uploaded_files"])
        shared_dict["times"]+=('\\Quality Filter('+str(round(qual_end-qual_start, 3)))
        shared_dict["times"]+=')'
        times = shared_dict["times"]
        return render_template("preview.html", times=times, num_rounds = shared_dict["num_rounds"], len_uploaded_files=len_uploaded_files, rounds=range(shared_dict["num_rounds"]), uploaded_files=shared_dict["uploaded_files"], all_high_qual_seqs=shared_dict["all_high_qual_seqs"], origonal_seqs=shared_dict["origonal_seqs"], quality_cutoff=shared_dict["quality_cutoff"])


@app.route('/preview', methods=["GET", "POST"])
def preview():
    if request.method == 'POST':

        preview_start = time.time()

        if not request.form.get("options"):
            return error_page("you need to enter an analysis type")
        if not request.form.get("rc"):
            return error_page("you need to enter if you want the reverse complement or not")

        trim_type = request.form.get("position_seq")
        shared_dict["trim_type"] = trim_type

        if trim_type == 'pos':

            if not request.form.get("start_stats"):
                return error_page("you need to enter a start position")
            if not request.form.get("end_stats"):
                return error_page("you need to enter an end position")

            pos_start = int(request.form.get("start_stats"))
            pos_end = int(request.form.get("end_stats"))
            if pos_start>=pos_end:
                return error_page("the start position must be less than the end position")
            if pos_start == 0 or pos_end == 0:
                return error_page("the start position and end position must be non zero")
            shared_dict["num_start"] = pos_start-1
            shared_dict["num_end"] = pos_end


        else:
            if trim_type == 'seq':

                if not request.form.get("start_seq"):
                    return error_page("you need to enter a start sequence")
                if not request.form.get("end_seq"):
                    return error_page("you need to enter an end sequence")



                seq_start = request.form.get("start_seq")
                seq_end = request.form.get("end_seq")

                allowed_trim = ['A', 'C', 'G', 'T']
                for letter in str(seq_start):
                    if letter not in allowed_trim:
                        return error_page("the start sequence contain only ACGT")

                for letter in str(seq_end):
                    if letter not in allowed_trim:
                        return error_page("the end sequence must contain only ACGT")

                shared_dict["seq_start"] = seq_start
                shared_dict["seq_end"] = seq_end

            else:
                return error_page("you need to select a type of trimmer")

        uploaded_files = []
        for file in glob.glob("*.fastq"):
            uploaded_files.append([get_uploaded_files(file)])

        shared_dict["option"] = request.form["options"]

        shared_dict["rc"] = request.form["rc"]

        if shared_dict["trim_type"] == 'pos':
            shared_dict["all_high_qual_seqs"], shared_dict["origonal_seqs"] = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["num_start"], shared_dict["num_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))
        if shared_dict["trim_type"] == 'seq':
            shared_dict["all_high_qual_seqs"], shared_dict["origonal_seqs"] = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["seq_start"], shared_dict["seq_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))


        if len(shared_dict["all_high_qual_seqs"][0])==0:
                        return error_page("please choose a lower quality threshold")
        min_len = len(shared_dict["all_high_qual_seqs"][0][0])

        for round in shared_dict["all_high_qual_seqs"]:
            for seq in round:
                if len(seq)<min_len:
                    min_len = len(seq)

        if (shared_dict["num_end"]-shared_dict["num_start"])>min_len:
            shared_dict["num_end"] = shared_dict["num_start"]+min_len

        if shared_dict["option"] == 'option1':
            preview_end = time.time()
            diff = preview_end-preview_start
            shared_dict["times"]+=('\\Preview('+str(diff)[:5])
            shared_dict["times"]+=')'
            times = shared_dict["times"]
            return render_template("data_analysis.html", times=times, num_rounds = shared_dict["num_rounds"], num_files = len(shared_dict["all_high_qual_seqs"]), num_start=shared_dict["num_start"], origonal_seqs=shared_dict["origonal_seqs"], rounds=range(shared_dict["num_rounds"]), num_end=shared_dict["num_end"], rc=shared_dict["rc"], all_high_qual_seqs=shared_dict["all_high_qual_seqs"], quality_cutoff=shared_dict["quality_cutoff"], uploaded_files=shared_dict["uploaded_files"])

        if shared_dict["option"] == 'option2':
            preview_end = time.time()
            diff = preview_end-preview_start
            shared_dict["times"]+=('\\Preview('+str(diff)[:5])
            shared_dict["times"]+=')'
            times = shared_dict["times"]
            return render_template("family_analysis.html", times=times, num_rounds = shared_dict["num_rounds"], num_files = len(shared_dict["all_high_qual_seqs"]),origonal_seqs=shared_dict["origonal_seqs"], num_start=shared_dict["num_start"], num_end=shared_dict["num_end"], rc=shared_dict["rc"], quality_cutoff=shared_dict["quality_cutoff"], uploaded_files=shared_dict["uploaded_files"], rounds = range(shared_dict["num_rounds"]))

    if request.method == "GET":
        times = shared_dict["times"]
        return render_template("preview.html", times=times, num_rounds = shared_dict["num_rounds"], len_uploaded_files=len_uploaded_files, rounds=range(shared_dict["num_rounds"]), uploaded_files=shared_dict["uploaded_files"], all_high_qual_seqs=shared_dict["all_high_qual_seqs"], origonal_seqs=shared_dict["origonal_seqs"], quality_cutoff=shared_dict["quality_cutoff"])



@app.route('/filter_qual', methods=["GET", "POST"])
def filter_qual():
    if request.method == 'POST':
        return render_template("quality.html", origonal_seqs=shared_dict["origonal_seqs"], rounds=range(len(shared_dict["all_high_qual_seqs"])), num_start=shared_dict["num_start"], num_end=shared_dict["num_end"], rc=shared_dict["rc"], all_high_qual_seqs=shared_dict["all_high_qual_seqs"], quality_cutoff=shared_dict["quality_cutoff"], uploaded_files=shared_dict["uploaded_files"])

    if request.method == "GET":
        return render_template("data_analysis.html", num_start=shared_dict["num_start"], origonal_seqs=shared_dict["origonal_seqs"], rounds=range(shared_dict["num_rounds"]), num_end=shared_dict["num_end"], rc=shared_dict["rc"], all_high_qual_seqs=shared_dict["all_high_qual_seqs"], quality_cutoff=shared_dict["quality_cutoff"], uploaded_files=shared_dict["uploaded_files"])

@app.route('/selection_statistics', methods=["GET", "POST"])
def selection_statistics():
    if request.method == 'POST':
        ss_start = time.time()
        global all_lens
        all_lens = []


        global all_seqs
        all_seqs = []

        for file in range(len(shared_dict["all_high_qual_seqs"])):
            global seqs
            seqs = shared_dict["all_high_qual_seqs"][file][1::4]
            qualities = shared_dict["all_high_qual_seqs"][file][3::4]

            acceptable_quals_1_percent_error = [5,6,7,8,9,':',';', '<', '=', '>', '?', '@', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K']

            all_lens.append([len(shared_dict["origonal_seqs"][file]), len(shared_dict["all_high_qual_seqs"][file]), len(list(set(shared_dict["all_high_qual_seqs"][file]))), round(len(list(set(shared_dict["all_high_qual_seqs"][file])))/len(shared_dict["all_high_qual_seqs"][file]), 3)])
        ss_end = time.time()
        shared_dict["times"]+=('\\Selection statistics('+str(round(ss_end-ss_start, 3)))
        shared_dict["times"]+=')'
        times = shared_dict["times"]
        return render_template("selection_statistics.html", times=times, all_lens=all_lens, rounds = range(len(all_lens)))

    if request.method == "GET":
        return render_template("data_analysis.html", num_start=shared_dict["num_start"], origonal_seqs=shared_dict["origonal_seqs"], rounds=range(shared_dict["num_rounds"]), num_rounds=shared_dict["num_rounds"], num_end=shared_dict["num_end"], rc=shared_dict["rc"], all_high_qual_seqs=shared_dict["all_high_qual_seqs"], quality_cutoff=shared_dict["quality_cutoff"], uploaded_files=shared_dict["uploaded_files"])


@app.route('/clustering', methods=["GET", "POST"])
def clustering():
    if request.method == 'POST':
        clustering_start = time.time()

        clustering_num = 0

        if shared_dict["top"] == []:
            return error_page("Please run Cluster Peak Finder first")

        labels = shared_dict["all_labels"][0]

        clustering_end = time.time()
        shared_dict["times"]+=('\\Clustering('+str(round(clustering_end-clustering_start, 3)))
        shared_dict["times"]+=')'
        times = shared_dict["times"]

        return render_template("cluster_seqs.html",times=times,num_clustering=clustering_num, num_rounds=shared_dict["num_rounds"], rounds=range(len(shared_dict["all_high_qual_seqs"])), all_labels=shared_dict["all_labels"],quality_cutoff=shared_dict["quality_cutoff"], all_high_qual_seqs=shared_dict["all_high_qual_seqs"], num_all_high_qual_seqs=len(shared_dict["all_high_qual_seqs"][0]))

    if request.method == "GET":
        return render_template("family_analysis.html")

@app.route('/most_abundant', methods=["GET", "POST"])
def most_abundant():
    if request.method == 'POST':
        ma_start = time.time()
        if not request.form.get("most_abundant"):
            return error_page("you need to enter a number of sequences")

        if request.form.get("most_abundant").isnumeric() == False:
            return error_page("you need to enter a number of sequences")

        num_analyze = int(request.form.get("most_abundant"))
        global analyze
        analyze = num_analyze

        if shared_dict["trim_type"] == 'pos':
            high_qual_seqs_not_unique, origonal_seqs_not_unique = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["num_start"], shared_dict["num_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))
        if shared_dict["trim_type"] == 'seq':
            high_qual_seqs_not_unique, origonal_seqs_not_unique = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["seq_start"], shared_dict["seq_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))

        counted_data = dict(collections.Counter(high_qual_seqs_not_unique[-1]))
        raw_data_list = [(k, v) for k, v in counted_data.items()]
        global top_N_sequences
        top_N_sequences = sorted(raw_data_list, key = lambda x: x[1], reverse=True)[:analyze]

        global most_abundant_seqs
        global most_abundant_nums
        most_abundant_seqs = []
        most_abundant_nums = []
        for seq in top_N_sequences:
            most_abundant_seqs.append([seq[0].strip('\n'), seq[1]])

        ma_end = time.time()
        shared_dict["times"]+=('\\Most Abundant('+str(round(ma_end-ma_start, 3)))
        shared_dict["times"]+=')'
        times = shared_dict["times"]

        shared_dict["most_abundant_seqs"] = most_abundant_seqs

        return render_template("most_abundant.html", times=times,top_N_sequences=sorted(raw_data_list, key = lambda x: x[1], reverse=True), quality_cutoff=shared_dict["quality_cutoff"], most_abundant_seqs=most_abundant_seqs, most_abundant_nums=most_abundant_nums,num_seqs=len(shared_dict["all_high_qual_seqs"][0]), analyze=analyze)

    if request.method == "GET":
        return render_template("data_analysis.html", num_start=shared_dict["num_start"])

@app.route('/download_abund', methods=["GET", "POST"])
def download_abund():
    if request.method == 'POST':
        da_start = time.time()
        input_data = shared_dict["most_abundant_seqs"]
        output_data = ''
        output_data+='Sequence,Reads\n'
        for line in range(min(1000, len(input_data))):
                output_data+=str(input_data[line][0]).strip('\n')+','+str(input_data[line][1]).strip('\n')+'\n'

        response = make_response(output_data)

        response.headers["Content-Disposition"] = "attachment; filename=most_abundant_sequences.csv"
        da_end = time.time()
        shared_dict["times"]+=('\\Download Abundant('+str(round(da_end-da_start, 3)))
        shared_dict["times"]+=')'
        times = shared_dict["times"]
        return response

    if request.method == "GET":
        return render_template("data_analysis.html", times=times, num_start=shared_dict["num_start"])


@app.route('/most_abundant_fams', methods=["GET", "POST"])
def most_abundant_fams():
    if request.method == 'POST':
        maf_start = time.time()
        if not request.form.get("abund_fams"):
            return error_page("you need to enter a max number of clusters")
        if not request.form.get("fam_abun"):
            return error_page("you need to enter a number of clusters")
        global num_clusters
        num_clusters = int(request.form.get("abund_fams"))
        shared_dict["num_clusters"] = num_clusters
        if num_clusters>50:
            return error_page("The max number of clusters is 50")
        global clusters_to_track
        clusters_to_track = int(request.form.get("fam_abun"))
        if num_clusters<1+clusters_to_track:
            clusters_to_track = num_clusters

        maf_start = time.time()

        counted_data = dict(collections.Counter(shared_dict["all_high_qual_seqs"][-1]))
        raw_data_list = [(k, v) for k, v in counted_data.items()]
        top_N_sequences = sorted(raw_data_list, key = lambda x: x[1], reverse=True)

        shared_dict["all_labels"] = get_clusters(num_clusters)
        labels = shared_dict["all_labels"][0]

        top = []
        used_clusters = []
        while len(top)<(clusters_to_track-1):
            for seq in range(len(top_N_sequences)):
                if labels[seq] not in used_clusters:
                    used_clusters.append(labels[seq])
                    top.append([top_N_sequences[seq], labels[seq]])

        shared_dict["top"] = top
        maf_end = time.time()
        shared_dict["times"]+=('\\Most Abundant Families('+str(round(maf_end-maf_start, 3)))
        shared_dict["times"]+=')'
        times = shared_dict["times"]

        shared_dict["top_clusters"]=top
        quality_cutoff=shared_dict["quality_cutoff"]

        return render_template("most_abundant_fams.html",times=times, num_clusters=num_clusters, top=top, quality_cutoff=shared_dict["quality_cutoff"], clusters_to_track=clusters_to_track)

@app.route('/dist_btwn', methods=["GET", "POST"])
def dist_btwn():
    if request.method == 'POST':
        db_start = time.time()
        if not request.form.get("dist_clusters"):
            return error_page("you need to enter a max number of clusters")
        if not request.form.get("dist_num"):
            return error_page("you need to enter a number of clusters")
        global num_clusters
        num_clusters = int(request.form.get("dist_clusters"))
        if num_clusters>50:
            return error_page("The max number of clusters is 50")

        global clusters_to_track
        clusters_to_track = int(request.form.get("dist_num"))
        if clusters_to_track<2:
            return error_page("You must enter 2 or more clusters to compare")
        if num_clusters<1+clusters_to_track:
            clusters_to_track = num_clusters

        shared_dict["all_labels"] = get_clusters(num_clusters)

        if shared_dict["trim_type"] == 'pos':
            high_qual_seqs_not_unique, origonal_seqs_not_unique = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["num_start"], shared_dict["num_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))
        if shared_dict["trim_type"] == 'seq':
            high_qual_seqs_not_unique, origonal_seqs_not_unique = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["seq_start"], shared_dict["seq_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))

        round = -1
        top_sequences = []
        for cluster in range(1+max(shared_dict["all_labels"][round])):
            cluster_seqs = []
            occurances = np.where(np.array(shared_dict["all_labels"][round]) == cluster)[0]
            max_counts = 0
            max_seq = ''
            for oc in occurances:
                counts = origonal_seqs_not_unique[round].count(str(shared_dict["all_high_qual_seqs"][round][oc]))
                if counts>max_counts:
                    max_counts = counts
                    max_seq = str(shared_dict["all_high_qual_seqs"][round][oc])
            top_sequences.append([max_seq, max_counts, cluster])

        top = sorted(top_sequences, key = lambda x: x[1], reverse=True)[:clusters_to_track]

        shared_dict["top_clusters"]=top

        all_diffs = []
        for seq1 in top:
            overlap=[]
            for seq2 in top:
                overlap.append((difflib.SequenceMatcher(None, seq1[0], seq2[0]).ratio())*len(seq1[0]))
            all_diffs.append(overlap)
        shared_dict["all_diffs"]=all_diffs

        db_end = time.time()
        shared_dict["times"]+=('\\Distance Between('+str(round(db_end-db_start, 3)))
        shared_dict["times"]+=')'
        times = shared_dict["times"]

        return render_template("dist_btwn.html", times=times, overlap=overlap, top_clusters=shared_dict["top_clusters"], all_diffs=all_diffs, top=shared_dict["top_clusters"], num_clusters=num_clusters, quality_cutoff=shared_dict["quality_cutoff"], clusters_to_track=clusters_to_track)

@app.route('/Cluster_Comparer.png')
def plot_btwn():
    fig = create_figure_btwn()
    output = io.BytesIO()
    FigureCanvas(fig).print_png(output)
    FigureCanvas(fig).print_png(output)
    return Response(output.getvalue(), mimetype='image/png')

def create_figure_btwn():
    fig = Figure()
    axis = fig.add_subplot(1, 1, 1)

    dis_matrix = pairwise_distances(np.array(shared_dict["all_diffs"]))
    mds_model = manifold.MDS(n_components = 2, random_state = 123, dissimilarity = 'precomputed')

    mds_fit = mds_model.fit(dis_matrix)
    mds_coords = mds_model.fit_transform(dis_matrix)

    labels = []
    for seq in range(len(dis_matrix)):
        labels.append('Seq'+str(seq))

    for el in range(len(mds_coords)):
        axis.scatter(mds_coords[el,0],mds_coords[el,1], s=40, marker = "o")

    n = list(range(1, 1+len(mds_coords)))
    for i, txt in enumerate(n):
        axis.annotate(txt, xy=(mds_coords[:,0][i], mds_coords[:,1][i]), xytext=(mds_coords[:,0][i]+3, mds_coords[:,1][i]+3))


    counter_60 = 0
    counter_70 = 0
    counter_80 = 0
    counter_90 = 0

    for el1 in range(len(mds_coords)):
        for el in range(len(mds_coords)):
                x_values = [mds_coords[el,0], mds_coords[el1,0]]
                y_values = [mds_coords[el,1], mds_coords[el1,1]]
                if shared_dict["all_diffs"][el][el1]/shared_dict["all_diffs"][el][0] <.4:
                    if counter_60 == 0:
                        axis.plot(x_values, y_values, color = 'red', label = '<40% similar')
                        counter_60+=1
                    else:
                        axis.plot(x_values, y_values, color = 'red')
                else:
                    if shared_dict["all_diffs"][el][el1]/shared_dict["all_diffs"][el][0] <.6:
                        if counter_70 ==0:
                            axis.plot(x_values, y_values, color = 'orange', label = '40-60% similar')
                            counter_70+=1
                        else:
                            axis.plot(x_values, y_values, color = 'orange')
                    else:
                        if shared_dict["all_diffs"][el][el1]/shared_dict["all_diffs"][el][0] <.8:
                            if counter_80 == 0:
                                axis.plot(x_values, y_values, color = 'green', label = '60-80% similar')
                                counter_80+=1
                            else:
                                axis.plot(x_values, y_values, color = 'green')
                        else:
                            if counter_90 == 0:
                                axis.plot(x_values, y_values, color = 'black', label = '>80% similar')
                                counter_90+=1
                            else:
                                axis.plot(x_values, y_values, color = 'black')

    axis.legend()

    axis.set_xlabel('First Dimension')
    axis.set_ylabel('Second Dimension')
    axis.set_title('Dissimilarity among most abundant clusters')

    axis.set_xlim(min(mds_coords[:,0])-5, 5+max(mds_coords[:,0]))
    axis.set_ylim(min(mds_coords[:,1])-5, 5+max(mds_coords[:,1]))

    return fig


@app.route('/track_seqs', methods=["GET", "POST"])
def track_seqs():

    if request.method == 'POST':
        if not request.form.get("track_seqs"):
            return error_page("you need to enter a number of sequences to track")
        if request.form.get("track_seqs").isnumeric() == False:
            return error_page("you need to enter a number of sequences to track")

        track_start = time.time()
        global analyze
        analyze = int(request.form.get("track_seqs"))


        if shared_dict["trim_type"] == 'pos':
            high_qual_seqs_not_unique, origonal_seqs_not_unique = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["num_start"], shared_dict["num_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))
        if shared_dict["trim_type"] == 'seq':
            high_qual_seqs_not_unique, origonal_seqs_not_unique = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["seq_start"], shared_dict["seq_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))


        counted_data = dict(collections.Counter(high_qual_seqs_not_unique[-1]))
        raw_data_list = [(k, v) for k, v in counted_data.items()]
        top_sequences = sorted(raw_data_list, key = lambda x: x[1], reverse=True)[:analyze]

        counting = []

        counts_from_prior_rounds = []
        for round in range(len(shared_dict["all_high_qual_seqs"])):

            counting.append(shared_dict["origonal_seqs"][round][0])

            counts = []
            for top in top_sequences:
                count = 0
                for seq in shared_dict["all_high_qual_seqs"][round]:
                    if top[0] in seq:
                        count+=1
                counts.append(count)
            counts_from_prior_rounds.append(counts)

        all_tracking_nums = []
        all_tracking_counts = []
        seq = 1
        for seq in range(analyze):
            tn_log = []
            tn_counts = []
            for round in range(len(counts_from_prior_rounds)):
                if counts_from_prior_rounds[round][seq]>0:
                    tn_log.append(counts_from_prior_rounds[round][seq])
                    tn_counts.append(counts_from_prior_rounds[round][seq])
                else:
                    tn_log.append(1)
                    tn_counts.append(counts_from_prior_rounds[round][seq])
            #all_tracking_nums.append(np.log(tn_log))
            all_tracking_nums.append(list(np.log(tn_log)))
            all_tracking_counts.append(np.array(tn_counts)/len(high_qual_seqs_not_unique[round]))

        track_end = time.time()
        diff = track_end-track_start
        shared_dict["times"]+=('\\Track Seqs('+str(diff)[:5])
        shared_dict["times"]+=')'
        times = shared_dict["times"]


        shared_dict["all_tracking_nums"]=all_tracking_nums
        shared_dict["all_tracking_counts"]=all_tracking_counts

        return render_template("track_seqs.html", times=times, fn=shared_dict["filenames"],counts_from_prior_rounds=counts_from_prior_rounds, ex=1, high_qual_seqs_not_unique=high_qual_seqs_not_unique, top_n=top_sequences, all_tracking_nums=all_tracking_nums, all_tracking_counts=all_tracking_counts, quality_cutoff=shared_dict["quality_cutoff"], num_track =analyze)

    if request.method == "GET":
        return render_template("data_analysis.html", num_start=shared_dict["num_start"])


@app.route('/Tracker_Sequences.png')
def plot_png_seqs():
    fig = create_figure()
    output = io.BytesIO()
    FigureCanvas(fig).print_png(output)
    FigureCanvas(fig).print_png(output)
    return Response(output.getvalue(), mimetype='image/png')

@app.route('/Tracker_Families.png')
def plot_png_fams():
    fig = create_figure()
    output = io.BytesIO()
    FigureCanvas(fig).print_png(output)
    FigureCanvas(fig).print_png(output)
    return Response(output.getvalue(), mimetype='image/png')

def create_figure():
    fig = Figure()
    axis = fig.add_subplot(1, 1, 1)

    for seq in range(len(shared_dict["all_tracking_counts"])):
        xs = range(1,1+len(shared_dict["all_tracking_counts"][seq]))
        ys = shared_dict["all_tracking_counts"][seq]
        axis.plot(xs, ys, marker = 'o', label = 'S'+str(1+seq))
        axis.xaxis.set_ticks(np.arange(1, 1+len(shared_dict["all_tracking_counts"][seq]), 1))
    axis.set_xlabel('Round')
    axis.set_ylabel('Abundance')
    fig.legend(bbox_to_anchor=(0.32, .95))
    return fig

@app.route('/Tracker_Bars.png')
def plot_nums():
    fig = create_figure_nums()
    output = io.BytesIO()
    FigureCanvas(fig).print_png(output)
    FigureCanvas(fig).print_png(output)
    return Response(output.getvalue(), mimetype='image/png')

def create_figure_nums():
    fig = Figure()
    axis = fig.add_subplot(111, projection='3d')
    cols = ['purple', 'indigo', 'violet', 'blue', 'skyblue', 'green', 'navyblue', 'yellow', 'orange', 'red','purple', 'indigo', 'violet', 'blue', 'skyblue', 'green', 'navyblue', 'yellow', 'orange', 'red']
    c=[]
    xpos = []
    for i in range(len(shared_dict["all_tracking_nums"])):
        xpos+=[i+1]*len(shared_dict["all_tracking_nums"][i])
    ypos = []
    for i in range(len(shared_dict["all_tracking_nums"])):
        ypos += list(range(len(shared_dict["all_tracking_nums"][i])))
        c+=[cols[i]]*len(shared_dict["all_tracking_nums"][i])

    zpos = [0]*len(xpos)
    dx = [.7]*(len(xpos))
    dy = [.3]*(len(xpos))
    dz = []
    for i in range(len(shared_dict["all_tracking_nums"])):
        dz+=shared_dict["all_tracking_nums"][i]

    axis.bar3d(xpos, ypos, zpos, dx, dy, dz, color=c)

    axis.set_xlabel('Sequence')
    axis.set_xticks(range(len(shared_dict["all_tracking_nums"])))
    axis.set_xticklabels(range(1, 1+len(shared_dict["all_tracking_nums"])))
    axis.set_ylabel('Round')
    axis.set_yticks(range(len(shared_dict["all_tracking_nums"][0])))
    axis.set_yticklabels(range(1, 1+len(shared_dict["all_tracking_nums"][0])))
    axis.set_zlabel('Log Counts')

    return fig

@app.route('/track_fams', methods=["GET", "POST"])
def track_fams():
    if request.method == 'POST':
        trackf_start = time.time()

        if shared_dict["top"] == []:
            return error_page("Please run Cluster Peak Finder first")

        if not request.form.get("to_track"):
            return error_page("you need to enter a number of families to track")

        to_track = int(request.form.get("to_track"))

        top = shared_dict["top"]
        quality_cutoff = shared_dict["quality_cutoff"]

        if shared_dict["trim_type"] == 'pos':
            high_qual_seqs_not_unique, origonal_seqs_not_unique = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["num_start"], shared_dict["num_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))
        if shared_dict["trim_type"] == 'seq':
            high_qual_seqs_not_unique, origonal_seqs_not_unique = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["seq_start"], shared_dict["seq_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))

        counting = []

        counts_from_prior_rounds = []
        for round in range(len(shared_dict["all_high_qual_seqs"])):

            counting.append(shared_dict["origonal_seqs"][round][0])

            counts = []
            for top_seq in top:
                count = 0
                for seq in shared_dict["all_high_qual_seqs"][round]:
                    if top_seq[0][0] in seq:
                        count+=1
                counts.append(count)
            counts_from_prior_rounds.append(counts)

        all_tracking_nums = []
        all_tracking_counts = []
        seq = 1
        for seq in range((to_track)):
            tn_log = []
            tn_counts = []
            for round in range(len(counts_from_prior_rounds)):
                if counts_from_prior_rounds[round][seq]>0:
                    tn_log.append(counts_from_prior_rounds[round][seq])
                    tn_counts.append(counts_from_prior_rounds[round][seq])
                else:
                    tn_log.append(1)
                    tn_counts.append(counts_from_prior_rounds[round][seq])
            #all_tracking_nums.append(np.log(tn_log))
            all_tracking_nums.append(list(np.log(tn_log)))
            all_tracking_counts.append(np.array(tn_counts)/len(high_qual_seqs_not_unique[round]))

        shared_dict["all_tracking_counts"] = []
        shared_dict["all_tracking_counts"] = all_tracking_counts

        shared_dict["all_tracking_nums"] = []
        shared_dict["all_tracking_nums"] = all_tracking_nums

        trackf_end = time.time()
        diff = trackf_end-trackf_start
        shared_dict["times"]+=('\\Track Families('+str(diff)[:5])
        shared_dict["times"]+=')'
        times = shared_dict["times"]

        return render_template("track_fams.html",times=times,all_tracking_nums=all_tracking_nums,top_clusters=top, all_tracking_counts=all_tracking_counts, quality_cutoff=quality_cutoff)

    if request.method == "GET":
        return render_template("family_analysis.html")

def create_figure_fams():
    fig = Figure()
    axis = fig.add_subplot(1, 1, 1)
    for seq in range(len(shared_dict["all_tracking_counts"])):
        xs = range(1,1+len(shared_dict["all_tracking_counts"][seq]))
        ys = shared_dict["all_tracking_counts"][seq]
        axis.plot(xs, ys, marker = 'o', label = 'Top_cluster_'+str(1+seq))
    axis.set_xlabel('Round')
    axis.set_ylabel('Abundance')
    fig.legend(bbox_to_anchor=(0.32, .95))
    return fig

@app.route('/plot_nums_fams.png')
def plot_nums_fams():
    fig = create_figure_nums_fams()
    output = io.BytesIO()
    FigureCanvas(fig).print_png(output)
    FigureCanvas(fig).print_png(output)
    return Response(output.getvalue(), mimetype='image/png')

def create_figure_nums_fams():
    fig = Figure()
    axis = fig.add_subplot(111, projection='3d')
    cols = ['purple', 'indigo', 'violet', 'blue', 'skyblue', 'green', 'navyblue', 'yellow', 'orange', 'red','purple', 'indigo', 'violet', 'blue', 'skyblue', 'green', 'navyblue', 'yellow', 'orange', 'red']
    c=[]
    xpos = []
    for i in range(len(shared_dict["all_tracking_nums"])):
        xpos+=[i+1]*len(shared_dict["all_tracking_nums"][i])
    ypos = []
    for i in range(len(shared_dict["all_tracking_nums"])):
        ypos += list(range(len(shared_dict["all_tracking_nums"][i])))
        c+=[cols[i]]*len(shared_dict["all_tracking_nums"][i])

    zpos = [0]*len(xpos)
    dx = [.7]*(len(xpos))
    dy = [.3]*(len(xpos))
    dz = []
    for i in range(len(shared_dict["all_tracking_nums"])):
        dz+=shared_dict["all_tracking_nums"][i]

    axis.bar3d(xpos, ypos, zpos, dx, dy, dz, color=c)

    axis.set_xlabel('Cluster')
    axis.set_xticks(range(len(shared_dict["all_tracking_nums"])))
    axis.set_xticklabels(range(0, 1+len(shared_dict["all_tracking_nums"])))
    axis.set_ylabel('Round')
    axis.set_yticks(range(len(shared_dict["all_tracking_nums"][0])))
    axis.set_yticklabels(range(1, 1+len(shared_dict["all_tracking_nums"][0])))
    axis.set_zlabel('Log Counts')

    return fig


@app.route('/seq_emergence', methods=["GET", "POST"])
def seq_emergence():
    if request.method == 'POST':
        if not request.form.get("seq_emergence"):
            return error_page("you need to enter a number of sequences")
        if not request.form.get("se_cmap"):
            return error_page("you need to enter a color scheme")

        global num_emerg
        num_emerg = int(request.form.get("seq_emergence"))

        global se_cmap
        se_cmap = request.form["se_cmap"]


        for file in range(shared_dict["num_rounds"]):

            top_n_seqs = []
            for round in range(len(shared_dict["all_high_qual_seqs"])):
                high_qual_seqs = shared_dict["all_high_qual_seqs"][round]

                all_counts = []
                for seq in high_qual_seqs:
                    all_counts.append(shared_dict["origonal_seqs"][round].count(seq))

                all_counts = np.array(all_counts)

                most_abundant_nums = all_counts.argsort()[-num_emerg:][::-1]
                global most_abundant_seqs
                most_abundant_seqs = []
                for ind in most_abundant_nums:
                    most_abundant_seqs.append([high_qual_seqs[ind], shared_dict["origonal_seqs"][round].count(high_qual_seqs[ind])])
                top_n_seqs.append(most_abundant_seqs)



            global total_heatmap_seqs
            total_heatmap_seqs = []

            for round2 in range(len(top_n_seqs)):
                for round1 in range(len(top_n_seqs)):
                    global all_overlap
                    all_overlap = []

                    for seq1 in range(len(top_n_seqs[round1])):
                        overlap = []
                        for seq2 in range(len(top_n_seqs[round2])):
                            overlap.append(difflib.SequenceMatcher(None, str(top_n_seqs[round1][seq1][0]),str(top_n_seqs[round2][seq2][0])).ratio())
                        all_overlap.append(overlap)
                    total_heatmap_seqs.append([all_overlap, round1, round2])


        return render_template("seq_emergence.html", se_cmap=se_cmap, quality_cutoff=shared_dict["quality_cutoff"], total_heatmap_seqs=total_heatmap_seqs, num_emerg =num_emerg, most_abundant_seqs=most_abundant_seqs)
    if request.method == "GET":
        return render_template("data_analysis.html", num_start=shared_dict["num_start"])

@app.route('/Sequence_Similarity_Finder.png')
def plot_se():
    fig = create_figure_se()
    output = io.BytesIO()
    FigureCanvas(fig).print_png(output)
    FigureCanvas(fig).print_png(output)
    return Response(output.getvalue(), mimetype='image/png')

def create_figure_se():
    fig, axs = plt.subplots(len(total_heatmap_seqs), figsize=(10, 5*len(total_heatmap_seqs)), facecolor='w', edgecolor='k')
    axs = axs.ravel()
    for i in range(len(total_heatmap_seqs)):
        if se_cmap == 'red':
            axs[i].imshow(np.array(total_heatmap_seqs[i][0]).transpose(), interpolation='nearest', cmap='YlOrRd')
        if se_cmap == 'blue':
            axs[i].imshow(np.array(total_heatmap_seqs[i][0]).transpose(), interpolation='nearest', cmap='YlGnBu')
        if se_cmap == 'bw':
            axs[i].imshow(np.array(total_heatmap_seqs[i][0]).transpose(), interpolation='nearest', cmap='Greys')


        axs[i].set_ylabel('Round ' + str(1+total_heatmap_seqs[i][1]))
        axs[i].set_xlabel('Round ' + str(1+total_heatmap_seqs[i][2]))
        for k in range(len(np.array(total_heatmap_seqs[i][0]).transpose())):
            for j in range(len(np.array(total_heatmap_seqs[i][0]).transpose())):
                text = axs[i].text(j, k, round(np.array(total_heatmap_seqs[i][0]).transpose()[k, j], 2),
                               ha="center", va="center", color="k")
    return fig



@app.route('/fam_emergence', methods=["GET", "POST"])
def fam_emergence():
    if request.method == 'POST':
        if not request.form.get("fam_emergence"):
            return error_page("you need to enter a number of clusters to track")
        if not request.form.get("track_fams"):
            return error_page("you need to enter a maximum number of clusters")

        global num_emerg
        num_emerg = int(request.form.get("fam_emergence"))
        global fams_emerg
        fams_emerg = int(request.form.get("track_fams"))

        if track_fams>50:
            return error_page("The max number of clusters is 50")

        if int(num_emerg)<int(fams_emerg):
            fams_emerg = num_emerg

        top_n_clusters = []

        for round_num in range(len(shared_dict["all_high_qual_seqs"])):
            high_quality_sequences = shared_dict["all_high_qual_seqs"][round_num]

            matrix = np.asarray([np.fromstring(str(s), dtype=np.uint8) for s in high_quality_sequences]);

            error = 10e10
            lowest_k = 10e10
            for k in range(1, 1+num_emerg):
                 print('Tried k='+str(k))
                 kmeans = KMeans(init="random", n_clusters=k,n_init=10,max_iter=300,random_state=42)
                 kmeans.fit(matrix)
                 kmeans.labels_[:]
                 if kmeans.inertia_<error:
                    error = kmeans.inertia_
                    lowest_k = k

            kmeans = KMeans(init="random", n_clusters=lowest_k,n_init=10,max_iter=300,random_state=42)
            kmeans.fit(matrix)

            global labels
            labels = kmeans.labels_[:]

            all_counts = []
            for seq in range(len(high_quality_sequences)):
                all_counts.append([high_quality_sequences[seq], shared_dict["origonal_seqs"][shared_dict["round_num"]].count(high_quality_sequences[seq]), labels[seq]])
            all_counts = sorted(all_counts, key=lambda x: x[0])


            used_clusters = []
            top_n = []
            seq = 1
            while len(used_clusters)<num_emerg:
                if all_counts[-seq][2] not in used_clusters:
                    used_clusters.append(all_counts[-seq][2])
                    top_n.append(all_counts[-seq])
                seq+=1
            top_n_clusters.append(top_n)

            global total_heatmap
            total_heatmap = []

            for round2 in range(len(top_n_clusters)):
                for round1 in range(len(top_n_clusters)):
                    global all_overlap
                    all_overlap = []

                    for seq1 in range(len(top_n_clusters[round1])):
                        overlap = []
                        for seq2 in range(len(top_n_clusters[round2])):
                            overlap.append(difflib.SequenceMatcher(None, str(top_n_clusters[round1][seq1][0]),str(top_n_clusters[round2][seq2][0])).ratio())
                        all_overlap.append(overlap)
                    total_heatmap.append([all_overlap, round1, round2])

        return render_template("fam_emergence.html", total_heatmap=total_heatmap, top_n_clusters=top_n_clusters, num_emerg =num_emerg, quality_cutoff=shared_dict["quality_cutoff"], fams_emerg=fams_emerg)

    if request.method == "GET":
        return render_template("family_analysis.html")

@app.route('/Cluster_Similarity_Finder.png')
def plot_fe():
    fig = create_figure_fe()
    output = io.BytesIO()
    FigureCanvas(fig).print_png(output)
    FigureCanvas(fig).print_png(output)
    return Response(output.getvalue(), mimetype='image/png')

def create_figure_fe():
    fig, axs = plt.subplots(len(total_heatmap), figsize=(10, 5*len(total_heatmap)), facecolor='w', edgecolor='k')
    axs = axs.ravel()
    for i in range(len(total_heatmap)):
        axs[i].imshow(np.array(total_heatmap[i][0]).transpose(), interpolation='nearest', cmap='YlOrRd')
        axs[i].set_ylabel('Round ' + str(1+total_heatmap[i][1]))
        axs[i].set_xlabel('Round ' + str(1+total_heatmap[i][2]))
        for k in range(len(np.array(total_heatmap[i][0]).transpose())):
            for j in range(len(np.array(total_heatmap[i][0]).transpose())):
                text = axs[i].text(j, k, round(np.array(total_heatmap[i][0]).transpose()[k, j], 2),
                               ha="center", va="center", color="k")

    return fig



@app.route('/conserved', methods=["GET", "POST"])
def conserved():
    if request.method == 'POST':
        con_start = time.time()
        if shared_dict["top"] == []:
            return error_page("Please run Cluster Peak Finder first")

        #if not request.form.get("conserved"):
        #    return error_page("you need to enter a conservation cutoff")
        if not request.form.get("hm_cmap"):
            return error_page("you need to enter a color scheme")
        if not request.form.get("rank"):
            return error_page("you need to enter a cluster rank")
        global hm_cmap
        hm_cmap = request.form["hm_cmap"]
        global rank
        rank = int(request.form.get("rank"))-1
        if rank>shared_dict["num_clusters"]:
            return error_page("Enter a lower rank")

        #in_conserved = request.form.get("conserved")
        #if in_conserved.isnumeric()==True:
        #    if float(in_conserved)>=1:
        #        global num_conserved
        #        num_conserved = int(request.form.get("conserved"))
        #else:
        #    return error_page("Please enter a conservation cutoff >1")

        top = shared_dict["top"]
        quality_cutoff = shared_dict["quality_cutoff"]

        cluster_to_use = top[rank][1]

        indices = [i for i, x in enumerate(shared_dict["all_labels"][-1]) if x == cluster_to_use]

        seqs_in_cluster = []
        for seq in indices:
            seqs_in_cluster.append(shared_dict["all_high_qual_seqs"][-1][seq])

        #conserved = []
        heatmap_data = []

        seqs_in_cluster = seqs_in_cluster[:100]

        for position in range(len(seqs_in_cluster[0])):
            options = []
            for seq in seqs_in_cluster:
                options.append(seq[position])

            raw_options = [(k, v) for k, v in collections.Counter(options).items()]
            nuc_list = []
            for element in raw_options:
                nuc_list.append(element[0])
            if 'A' not in nuc_list:
                raw_options.append(('A', 0))
            if 'C' not in nuc_list:
                raw_options.append(('C', 0))
            if 'G' not in nuc_list:
                raw_options.append(('G', 0))
            if 'T' not in nuc_list:
                raw_options.append(('T', 0))
            raw_options = sorted(raw_options, key=lambda x: x[0])

            #for element in raw_options:
            #    if num_conserved>1:
            #        if element[1]/len(seqs_in_cluster)>(num_conserved/100):
                        #conserved.append(raw_options)
            #    else:
            #        if element[1]/len(seqs_in_cluster)>(num_conserved):
            #            conserved.append(raw_options)

            heatmap_data_position = []
            for element in raw_options:
                heatmap_data_position.append(element[1]/len(seqs_in_cluster))
            heatmap_data.append(heatmap_data_position)

        heatmap_data = heatmap_data[:-1]

        counter = -1
        for element in heatmap_data:
            counter+=1
            if max(element) == 1:
                heatmap_data[counter] = [0,0,0,0]

        con_end = time.time()
        diff = con_end-con_start
        shared_dict["times"]+=('\\Conserved('+str(diff)[:5])
        shared_dict["times"]+=')'
        times = shared_dict["times"]

        shared_dict["heatmap_data"] = heatmap_data

        lst = heatmap_data
        dfObj = pd.DataFrame(lst, columns = ['A' , 'C', 'G' , 'T'], index=range(len(lst)))
        shared_dict["dfObj"] = dfObj

        return render_template("conserved_positions.html", all_scores = shared_dict["all_scores"], times=times, heatmap_data=heatmap_data, raw_options=raw_options, hm_cmap=hm_cmap, rank=rank)
        #return render_template("conserved_positions.html", times=times, heatmap_data=heatmap_data, conserved=conserved, raw_options=raw_options, hm_cmap=hm_cmap, rank=rank, num_conserved=num_conserved)

    if request.method == "GET":
        return render_template("family_analysis.html")



@app.route('/Conservation_Finder.png')
def plot_hm():
    fig = create_figure_hm()
    output = io.BytesIO()
    FigureCanvas(fig).print_png(output)
    FigureCanvas(fig).print_png(output)
    return Response(output.getvalue(), mimetype='image/png')

@app.route('/Conservation_Logo.png')
def plot_cl():
    fig = create_figure_logo()
    output = io.BytesIO()
    FigureCanvas(fig).print_png(output)
    return Response(output.getvalue(), mimetype='image/png')

class Scale(matplotlib.patheffects.RendererBase):
    def __init__(self, sx, sy=None):
        self._sx = sx
        self._sy = sy

    def draw_path(self, renderer, gc, tpath, affine, rgbFace):
        affine = affine.identity().scale(self._sx, self._sy)+affine
        renderer.draw_path(gc, tpath, affine, rgbFace)

def create_figure_logo():
    all_scores = []
    nts = ['A','C','G','T']
    for pos in shared_dict["heatmap_data"]:
        score = []
        counter = 0
        for nt in pos:
            score.append((nts[counter], nt))
            counter+=1
        all_scores.append(score)

    shared_dict["all_scores"] = all_scores
    COLOR_SCHEME = {'G': 'orange',
                'A': 'red',
                'C': 'blue',
                'T': 'darkgreen'}
    BASES = list(COLOR_SCHEME.keys())
    fig = plt.figure()
    fig.set_size_inches(len(all_scores),2.5)
    ax = fig.add_subplot(111)
    ax.set_xticks(range(len(all_scores)))

    xshift = 0
    trans_offset = transforms.offset_copy(ax.transAxes,
                                      fig=fig,
                                      x=0,
                                      y=0,
                                      units='points')


    for scores in all_scores:
        yshift = 0
        for base, score in scores:
            txt = ax.text(0,
                          0,
                          base,
                          transform=trans_offset,
                          fontsize=80,
                          color=COLOR_SCHEME[base],
                          weight='bold',
                          ha='center',
                          family='sans-serif'
                          )
            txt.set_clip_on(False)
            txt.set_path_effects([Scale(1.0, score)])
            fig.canvas.draw()
            window_ext = txt.get_window_extent(txt._renderer)
            yshift = window_ext.height*score
            trans_offset = transforms.offset_copy(txt._transform, fig=fig, y=yshift, units='points')
            #trans_offset = transforms.offset_copy(txt._transform, fig=fig, y=y, units='points')
        xshift += window_ext.width
        trans_offset = transforms.offset_copy(ax.transAxes, fig=fig, x=xshift, units='points')


    ax.set_yticks(range(0,3))


    seaborn.despine(ax=ax, offset=30, trim=True)
    ax.set_xticklabels(range(1,len(all_scores)+1), rotation=90)
    ax.set_yticklabels(np.arange(0,3,1))
    return fig

def create_figure_hm():
    fig, ax = plt.subplots(figsize = (12,6))
    if hm_cmap == 'red':
        im = ax.imshow(np.array(shared_dict["heatmap_data"]).transpose(), interpolation='nearest', cmap='YlOrRd')
        y_label_list = ['', 'A', 'C', 'G', 'T', '']
        ax.set_yticks(range(-1, 1+len(shared_dict["heatmap_data"][0])))
        ax.set_yticklabels(y_label_list)
        ax.set_xticks(range(len(shared_dict["heatmap_data"])))
        ax.set_xticklabels(range(1, len(shared_dict["heatmap_data"])+1))

    if hm_cmap == 'blue':
        im = ax.imshow(np.array(shared_dict["heatmap_data"]).transpose(), interpolation='nearest', cmap='YlGnBu')
        y_label_list = ['', 'A', 'C', 'G', 'T', '']
        ax.set_yticks(range(-1, 1+len(shared_dict["heatmap_data"][0])))
        ax.set_yticklabels(y_label_list)
        ax.set_xticks(range(len(shared_dict["heatmap_data"])))
        ax.set_xticklabels(range(1, len(shared_dict["heatmap_data"])+1))


    if hm_cmap == 'bw':
        im = ax.imshow(np.array(shared_dict["heatmap_data"]).transpose(), interpolation='nearest', cmap='Greys')
        y_label_list = ['', 'A', 'C', 'G', 'T', '']
        ax.set_yticks(range(-1, 1+len(shared_dict["heatmap_data"][0])))
        ax.set_yticklabels(y_label_list)
        ax.set_xticks(range(len(shared_dict["heatmap_data"])))
        ax.set_xticklabels(range(1, len(shared_dict["heatmap_data"])+1))

    cbar = fig.colorbar(im, orientation='horizontal')
    #cbar.set_clim(0, 1)
    return fig

@app.route("/structure", methods=["GET", "POST"])
def structure():
    if request.method == "POST":
        if not request.form.get("ss_emergence"):
            return error_page("you need to enter a structure")
        global structure
        structure =  str(request.form.get("ss_emergence"))
        allowed = ['.', '(', ')']
        for letter in structure:
            if letter not in allowed:
                return error_page("the structure must contain only .,(,)")

        round = -1
        high_qual_seqs = shared_dict["all_high_qual_seqs"][round]
        all_counts = []
        for seq in high_qual_seqs:
            all_counts.append(shared_dict["origonal_seqs"][round].count(seq))

        all_counts = np.array(all_counts)

        most_abundant_nums = all_counts.argsort()[-2:][::-1]
        global most_abundant_seqs
        most_abundant_seqs = []
        for ind in most_abundant_nums:
            most_abundant_seqs.append(high_qual_seqs[ind])

        global consistant_ss
        consistant_ss = []
        global all_ss
        all_ss = []
        global seq_use
        seq_use = most_abundant_seqs[0]


        options = webdriver.ChromeOptions()
        options.add_argument('--ignore-certificate-errors')
        options.add_argument("--test-type")
        options.binary_location = "/usr/bin/chromium"
        driver = webdriver.Chrome(chrome_options=options)

        driver.get("http://www.google.com")
        display.stop()


        #url = "http://rna.tbi.univie.ac.at//cgi-bin/RNAWebSuite/RNAfold.cgi"
        #br = mechanize.Browser()
        #br.set_handle_robots(False) # ignore robots
        #br.addheaders = [('User-agent', 'Mozilla/5.0 (X11; U; Linux i686; en-US; rv:1.9.0.1) Gecko/2008071615 Fedora/3.0.1-1.fc9 Firefox/3.0.1')]

        #br.open(url)
        #br.select_form(name="form")
        #br["SCREEN"] = seq_use
        #res = br.submit()
        #content = str(res.read())

  #energy = int(re.sub("[^0-9]", "", str(content).split('thermodynamic ensemble is')[1].split('kcal/mol')[0]))
  #structure = str(content).split('STRUCTURE=')[1].split('target')[0].split('"')[0]

        return render_template("structure.html", seq_use=seq_use, content=content, all_ss=all_ss, top_seqs=most_abundant_seqs, consistant_ss=consistant_ss, structure=structure, quality_cutoff=shared_dict["quality_cutoff"])

    if request.method == "GET":
        return render_template("data_analysis.html", structure=structure, quality_cutoff=shared_dict["quality_cutoff"], num_start=shared_dict["num_start"])


@app.route('/paths', methods=["GET", "POST"])
def paths():
    if request.method == 'POST':
        paths_start = time.time()
        if not request.form.get("seq1"):
            return error_page("you need to enter a seq1")
        if not request.form.get("seq2"):
            return error_page("you need to enter a seq2")

        global seq1
        seq1 = str(request.form.get("seq1").replace(' ', ''))
        global seq2
        seq2 = str(request.form.get("seq2").replace(' ', ''))

        total_counts_1 = 0
        for round in shared_dict["all_high_qual_seqs"]:
            total_counts_1+=round.count(seq1)
        total_counts_2 = 0
        for round in shared_dict["all_high_qual_seqs"]:
            total_counts_2+=round.count(seq2)
        if total_counts_1 == 0:
            return error_page("you need to enter a starting sequence in the filtered sequences - see Download Preprocessed data spreadsheet")
        if total_counts_2 == 0:
            return error_page("you need to enter an ending sequence in the filtered sequences - see Download Preprocessed data spreadsheet")

        same_positions = []
        for position in range(len(seq1)):
            if seq1[position]==seq2[position]:
                same_positions.append(position)

        intermediate = []

        for test_seq in list(set(shared_dict["all_high_qual_seqs"][-1])):
            count = 0
            for position in same_positions:
                if test_seq[position]==seq1[position]:
                    count+=1
            if count > len(same_positions)-1:
                intermediate.append([test_seq, shared_dict["origonal_seqs"][-1].count(test_seq)])
        paths_end = time.time()
        diff = paths_end-paths_start
        shared_dict["times"]+=('\\Intermediates('+str(diff)[:5])
        shared_dict["times"]+=')'
        times = shared_dict["times"]

        return render_template("paths.html", times=times, intermediate=intermediate, seq1=seq1, seq2=seq2)

    if request.method == "GET":
        return render_template("data_analysis.html", num_start=shared_dict["num_start"])


@app.route('/region_emergence', methods=["GET", "POST"])
def region_emergence():
    if request.method == 'POST':
        region_start = time.time()
        if not request.form.get("region_emergence"):
            return error_page("you need to enter a region")

        global region
        region = str(request.form.get("region_emergence"))
        allowed = ['A', 'C', 'G', 'T', 'R', 'Y', 'N']
        for letter in region:
            if letter not in allowed:
                return error_page("the region must contain only ACGT")

        global chosen_region
        chosen_region = 0

        R_regions=[]

        if 'R' in region:
            R_regions.append(region.replace('R', 'G'))
            R_regions.append(region.replace('R', 'A'))
        else:
            R_regions.append(region)

        Y_regions=[]


        for element in R_regions:
            if 'Y' in element:
                Y_regions.append(element.replace('Y', 'C'))
                Y_regions.append(element.replace('Y', 'T'))
            else:
                for element in R_regions:
                    Y_regions.append(element)

        N_regions = []

        for element in Y_regions:
            if 'N' in element:
                N_regions.append(element.replace('N', 'C'))
                N_regions.append(element.replace('N', 'T'))
                N_regions.append(element.replace('N', 'G'))
                N_regions.append(element.replace('N', 'A'))
            else:
                for element in Y_regions:
                    N_regions.append(element)

        global all_regions
        all_regions = list(set(N_regions))
        for element in all_regions:
            if 'N' in element or 'Y' in element or 'R' in element:
                all_regions.remove(elemenet)

        global all_regions_all_counts_rounds
        all_regions_all_counts_rounds = []

        global all_regions_all_counts_rounds_seqs
        all_regions_all_counts_rounds_seqs = []

        for reg in all_regions:


            for file in range(shared_dict["num_rounds"]):
                global all_seqs_with_region
                all_seqs_with_region = []
                global all_counts_rounds
                all_counts_rounds = []
                for round_num in range(len(shared_dict["all_high_qual_seqs"])):
                    counts_round = 0
                    total_c = 0
                    seqs_with_region = []
                    for seq in shared_dict["all_high_qual_seqs"][round_num]:
                        total_c+=1
                        if str(reg) in str(seq):
                            counts_round+=1
                            seqs_with_region.append(seq)
                    all_counts_rounds.append(counts_round/total_c)
                    all_seqs_with_region.append(seqs_with_region)
            all_regions_all_counts_rounds_seqs.append(all_seqs_with_region)
            all_regions_all_counts_rounds.append(all_counts_rounds)
        shared_dict["all_regions_all_counts_rounds_seqs"] = all_regions_all_counts_rounds_seqs
        region_end = time.time()
        shared_dict["times"]+=('\\Motif('+str(round(region_end-region_start, 3)))
        shared_dict["times"]+=')'
        times = shared_dict["times"]
        all_regions_all_seqs_with_region = shared_dict["all_regions_all_counts_rounds_seqs"]
        return render_template("region_emergence.html", all_regions_all_counts_rounds_seqs=all_regions_all_counts_rounds_seqs, times=times,num_all_regions = range(len(all_regions)), all_regions=all_regions, all_regions_all_counts_rounds=all_regions_all_counts_rounds, all_regions_all_seqs_with_region=all_regions_all_seqs_with_region, rounds = range(len(shared_dict["all_high_qual_seqs"])), all_seqs_with_region=all_seqs_with_region, quality_cutoff=shared_dict["quality_cutoff"], all_counts_rounds=all_counts_rounds, num_clustering =clustering, num_seqs_to_cluster = len(shared_dict["all_high_qual_seqs"]), region=region)
    if request.method == "GET":
        return render_template("data_analysis.html", num_start=shared_dict["num_start"])

def create_figure_region():
    fig = Figure()

    axis = fig.add_subplot(1, 1, 1)

    for element in range(len(all_regions_all_counts_rounds)):

        xs = range(1,1+len(all_regions_all_counts_rounds[element]))
        ys = all_regions_all_counts_rounds[element]
        axis.plot(xs, ys, marker = 'o', label = str(all_regions[element]))
        axis.xaxis.set_ticks(np.arange(1, 1+len(all_regions_all_counts_rounds[element]), 1))

    axis.set_xlabel('Round')
    axis.set_ylabel('Abundance')
    fig.legend(bbox_to_anchor=(0.32, .95))

    return fig

@app.route('/Sequence_Motif_Finder.png')
def plot_region():
    fig = create_figure_region()
    output = io.BytesIO()
    FigureCanvas(fig).print_png(output)
    fig = create_figure_region()
    output = io.BytesIO()
    FigureCanvas(fig).print_png(output)
    return Response(output.getvalue(), mimetype='image/png')

@app.route("/", methods=["GET", "POST"])
def home():
    if request.method == "GET":
        return redirect("/home")

@app.route("/home", methods=["GET", "POST"])
def home_route():
    if request.method == "GET":
        return render_template("home.html")

@app.route("/about", methods=["GET"])
def about():
    if request.method == "GET":
        return render_template("about.html")

@app.route("/start", methods=["POST"])
def start():
    if request.method == "POST":
        return redirect("/choose_rounds")

@app.route("/data_analysis", methods=["POST"])
def data_analysis():
    if request.method == "POST":
        return render_template("data_analysis.html", num_rounds = shared_dict["num_rounds"])

@app.route("/family_analysis", methods=["POST"])
def family_analysis():
    if request.method == "POST":
        return render_template("family_analysis.html", num_rounds = shared_dict["num_rounds"])

@app.route("/preview_back", methods=["POST"])
def preview_back():
    if request.method == "POST":
        return render_template("preview.html", num_start=shared_dict["num_start"], origonal_seqs=shared_dict["origonal_seqs"], rounds=range(shared_dict["num_rounds"]), num_rounds=shared_dict["num_rounds"], num_end=shared_dict["num_end"], rc=shared_dict["rc"], all_high_qual_seqs=shared_dict["all_high_qual_seqs"], quality_cutoff=shared_dict["quality_cutoff"], uploaded_files=shared_dict["uploaded_files"])


@app.route("/contact", methods=["GET"])
def contact():
    if request.method == "GET":
        return render_template("contact.html")

@app.route("/tutorial", methods=["GET"])
def tutorial():
    if request.method == "GET":
        return render_template("tutorial.html")

@app.route("/help", methods=["GET"])
def help():
    if request.method == "GET":
        return render_template("help.html")

@app.route("/resources", methods=["GET"])
def resources():
    if request.method == "GET":
        return render_template("resources.html")

@app.route('/pdf/<path:filename>', methods=['GET', 'POST'])
def download(filename):
    return send_from_directory(directory='pdf', filename=filename)


##----------------- CONSOLIDATE W LOOP------------------##
@app.route('/download_qual_pre_0', methods=["GET", "POST"])
def download_qual_pre_0():
    if request.method == 'POST':

        if shared_dict["trim_type"] == 'pos':
            high_qual_seqs_not_unique, origonal_seqs_not_unique = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["num_start"], shared_dict["num_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))
        if shared_dict["trim_type"] == 'seq':
            high_qual_seqs_not_unique, origonal_seqs_not_unique = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["seq_start"], shared_dict["seq_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))

        rounds = range(len(shared_dict["all_high_qual_seqs"]))

        input_data = high_qual_seqs_not_unique[0]
        output_data = ''
        output_data+='Sequence,Reads\n'

        raw_data_list = [(k, v) for k, v in collections.Counter(input_data).items()]
        top_sequences = sorted(raw_data_list, key = lambda x: x[1], reverse=True)

        for line in top_sequences:
            spaced = ''
            for l in range(len(line[0])):
                if (l+1)%3==0:
                    spaced+=str(line[0][l])+' '
                else:
                    spaced+=str(line[0][l])
            #output_data+=str(spaced).strip('\n')+','+str(line[1]).strip('\n')+'\n'
            output_data+=str(line[0]).strip('\n')+','+str(line[1]).strip('\n')+'\n'

        response = make_response(output_data)

        response.headers["Content-Disposition"] = "attachment; filename=High_Quality_Sequences_Round1.csv"
        return response

    if request.method == "GET":
        return render_template("quality.html")

@app.route('/download_qual_pre_1', methods=["GET", "POST"])
def download_qual_pre_1():
    if request.method == 'POST':

        if shared_dict["trim_type"] == 'pos':
            high_qual_seqs_not_unique, origonal_seqs_not_unique = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["num_start"], shared_dict["num_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))
        if shared_dict["trim_type"] == 'seq':
            high_qual_seqs_not_unique, origonal_seqs_not_unique = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["seq_start"], shared_dict["seq_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))

        rounds = range(len(shared_dict["all_high_qual_seqs"]))

        input_data = high_qual_seqs_not_unique[1]
        output_data = ''
        output_data+='Sequence,Reads\n'

        raw_data_list = [(k, v) for k, v in collections.Counter(input_data).items()]
        top_sequences = sorted(raw_data_list, key = lambda x: x[1], reverse=True)[:10000]

        for line in top_sequences:
            spaced = ''
            for l in range(len(line[0])):
                if (l+1)%3==0:
                    spaced+=str(line[0][l])+' '
                else:
                    spaced+=str(line[0][l])
            #output_data+=str(spaced.strip('\n'))+','+str(line[1])+'\n'
            output_data+=str(line[0]).strip('\n')+','+str(line[1]).strip('\n')+'\n'

        response = make_response(output_data)

        response.headers["Content-Disposition"] = "attachment; filename=High_Quality_Sequences_Round2.csv"
        return response

    if request.method == "GET":
        return render_template("quality.html")
@app.route('/download_qual_pre_2', methods=["GET", "POST"])
def download_qual_pre_2():
    if request.method == 'POST':

        if shared_dict["trim_type"] == 'pos':
            high_qual_seqs_not_unique, origonal_seqs_not_unique = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["num_start"], shared_dict["num_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))
        if shared_dict["trim_type"] == 'seq':
            high_qual_seqs_not_unique, origonal_seqs_not_unique = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["seq_start"], shared_dict["seq_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))

        rounds = range(len(shared_dict["all_high_qual_seqs"]))

        input_data = high_qual_seqs_not_unique[2]
        output_data = ''
        output_data+='Sequence,Reads\n'

        raw_data_list = [(k, v) for k, v in collections.Counter(input_data).items()]
        top_sequences = sorted(raw_data_list, key = lambda x: x[1], reverse=True)[:10000]

        for line in top_sequences:
            spaced = ''
            for l in range(len(line[0])):
                if (l+1)%3==0:
                    spaced+=str(line[0][l])+' '
                else:
                    spaced+=str(line[0][l])
            #output_data+=str(spaced.strip('\n'))+','+str(line[1])+'\n'
            output_data+=str(line[0]).strip('\n')+','+str(line[1]).strip('\n')+'\n'

        response = make_response(output_data)

        response.headers["Content-Disposition"] = "attachment; filename=High_Quality_Sequences_Round3.csv"
        return response

    if request.method == "GET":
        return render_template("quality.html")
@app.route('/download_qual_pre_3', methods=["GET", "POST"])
def download_qual_pre_3():
    if request.method == 'POST':

        if shared_dict["trim_type"] == 'pos':
            high_qual_seqs_not_unique, origonal_seqs_not_unique = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["num_start"], shared_dict["num_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))
        if shared_dict["trim_type"] == 'seq':
            high_qual_seqs_not_unique, origonal_seqs_not_unique = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["seq_start"], shared_dict["seq_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))

        rounds = range(len(shared_dict["all_high_qual_seqs"]))

        input_data = high_qual_seqs_not_unique[3]
        output_data = ''
        output_data+='Sequence,Reads\n'

        raw_data_list = [(k, v) for k, v in collections.Counter(input_data).items()]
        top_sequences = sorted(raw_data_list, key = lambda x: x[1], reverse=True)[:10000]

        for line in top_sequences:
            spaced = ''
            for l in range(len(line[0])):
                if (l+1)%3==0:
                    spaced+=str(line[0][l])+' '
                else:
                    spaced+=str(line[0][l])
            #output_data+=str(spaced.strip('\n'))+','+str(line[1])+'\n'
            output_data+=str(line[0]).strip('\n')+','+str(line[1]).strip('\n')+'\n'

        response = make_response(output_data)

        response.headers["Content-Disposition"] = "attachment; filename=High_Quality_Sequences_Round4.csv"
        return response

    if request.method == "GET":
        return render_template("quality.html")
@app.route('/download_qual_pre_4', methods=["GET", "POST"])
def download_qual_pre_4():
    if request.method == 'POST':

        if shared_dict["trim_type"] == 'pos':
            high_qual_seqs_not_unique, origonal_seqs_not_unique = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["num_start"], shared_dict["num_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))
        if shared_dict["trim_type"] == 'seq':
            high_qual_seqs_not_unique, origonal_seqs_not_unique = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["seq_start"], shared_dict["seq_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))

        rounds = range(len(shared_dict["all_high_qual_seqs"]))

        input_data = high_qual_seqs_not_unique[4]
        output_data = ''
        output_data+='Sequence,Reads\n'

        raw_data_list = [(k, v) for k, v in collections.Counter(input_data).items()]
        top_sequences = sorted(raw_data_list, key = lambda x: x[1], reverse=True)[:10000]

        for line in top_sequences:
            spaced = ''
            for l in range(len(line[0])):
                if (l+1)%3==0:
                    spaced+=str(line[0][l])+' '
                else:
                    spaced+=str(line[0][l])
            #output_data+=str(spaced.strip('\n'))+','+str(line[1])+'\n'
            output_data+=str(line[0]).strip('\n')+','+str(line[1]).strip('\n')+'\n'

        response = make_response(output_data)

        response.headers["Content-Disposition"] = "attachment; filename=High_Quality_Sequences_Round5.csv"
        return response

    if request.method == "GET":
        return render_template("quality.html")
@app.route('/download_qual_pre_5', methods=["GET", "POST"])
def download_qual_pre_5():
    if request.method == 'POST':

        rounds = range(len(shared_dict["all_high_qual_seqs"]))

        if shared_dict["trim_type"] == 'pos':
            high_qual_seqs_not_unique, origonal_seqs_not_unique = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["num_start"], shared_dict["num_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))
        if shared_dict["trim_type"] == 'seq':
            high_qual_seqs_not_unique, origonal_seqs_not_unique = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["seq_start"], shared_dict["seq_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))

        rounds = range(len(shared_dict["all_high_qual_seqs"]))

        input_data = high_qual_seqs_not_unique[5]
        output_data = ''
        output_data+='Sequence,Reads\n'

        raw_data_list = [(k, v) for k, v in collections.Counter(input_data).items()]
        top_sequences = sorted(raw_data_list, key = lambda x: x[1], reverse=True)[:10000]

        for line in top_sequences:
            spaced = ''
            for l in range(len(line[0])):
                if (l+1)%3==0:
                    spaced+=str(line[0][l])+' '
                else:
                    spaced+=str(line[0][l])
            output_data+=str(spaced.strip('\n'))+','+str(line[1])+'\n'
            #output_data+=str(line[0]).strip('\n')+','+str(line[1]).strip('\n')+'\n'

        response = make_response(output_data)

        response.headers["Content-Disposition"] = "attachment; filename=High_Quality_Sequences_Round6.csv"
        return response

    if request.method == "GET":
        return render_template("quality.html")
@app.route('/download_qual_pre_6', methods=["GET", "POST"])
def download_qual_pre_6():
    if request.method == 'POST':

        if shared_dict["trim_type"] == 'pos':
            high_qual_seqs_not_unique, origonal_seqs_not_unique = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["num_start"], shared_dict["num_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))
        if shared_dict["trim_type"] == 'seq':
            high_qual_seqs_not_unique, origonal_seqs_not_unique = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["seq_start"], shared_dict["seq_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))

        rounds = range(len(shared_dict["all_high_qual_seqs"]))

        input_data = high_qual_seqs_not_unique[6]
        output_data = ''
        output_data+='Sequence,Reads\n'

        raw_data_list = [(k, v) for k, v in collections.Counter(input_data).items()]
        top_sequences = sorted(raw_data_list, key = lambda x: x[1], reverse=True)[:10000]

        for line in top_sequences:
            spaced = ''
            for l in range(len(line[0])):
                if (l+1)%3==0:
                    spaced+=str(line[0][l])+' '
                else:
                    spaced+=str(line[0][l])
            #output_data+=str(spaced.strip('\n'))+','+str(line[1])+'\n'
            output_data+=str(line[0]).strip('\n')+','+str(line[1]).strip('\n')+'\n'

        response = make_response(output_data)

        response.headers["Content-Disposition"] = "attachment; filename=High_Quality_Sequences_Round7.csv"
        return response

    if request.method == "GET":
        return render_template("quality.html")
@app.route('/download_qual_pre_7', methods=["GET", "POST"])
def download_qual_pre_7():
    if request.method == 'POST':

        rounds = range(len(shared_dict["all_high_qual_seqs"]))

        if shared_dict["trim_type"] == 'pos':
            high_qual_seqs_not_unique, origonal_seqs_not_unique = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["num_start"], shared_dict["num_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))
        if shared_dict["trim_type"] == 'seq':
            high_qual_seqs_not_unique, origonal_seqs_not_unique = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["seq_start"], shared_dict["seq_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))

        rounds = range(len(shared_dict["all_high_qual_seqs"]))

        input_data = high_qual_seqs_not_unique[7]
        output_data = ''
        output_data+='Sequence,Reads\n'

        raw_data_list = [(k, v) for k, v in collections.Counter(input_data).items()]
        top_sequences = sorted(raw_data_list, key = lambda x: x[1], reverse=True)[:10000]

        for line in top_sequences:
            spaced = ''
            for l in range(len(line[0])):
                if (l+1)%3==0:
                    spaced+=str(line[0][l])+' '
                else:
                    spaced+=str(line[0][l])
            output_data+=str(spaced.strip('\n'))+','+str(line[1])+'\n'
            #output_data+=str(line[0]).strip('\n')+','+str(line[1]).strip('\n')+'\n'

        response = make_response(output_data)

        response.headers["Content-Disposition"] = "attachment; filename=High_Quality_Sequences_Round8.csv"
        return response

    if request.method == "GET":
        return render_template("quality.html")
@app.route('/download_qual_pre_8', methods=["GET", "POST"])
def download_qual_pre_8():
    if request.method == 'POST':

        if shared_dict["trim_type"] == 'pos':
            high_qual_seqs_not_unique, origonal_seqs_not_unique = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["num_start"], shared_dict["num_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))
        if shared_dict["trim_type"] == 'seq':
            high_qual_seqs_not_unique, origonal_seqs_not_unique = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["seq_start"], shared_dict["seq_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))

        rounds = range(len(shared_dict["all_high_qual_seqs"]))

        input_data = high_qual_seqs_not_unique[8]
        output_data = ''
        output_data+='Sequence,Reads\n'

        raw_data_list = [(k, v) for k, v in collections.Counter(input_data).items()]
        top_sequences = sorted(raw_data_list, key = lambda x: x[1], reverse=True)[:10000]

        for line in top_sequences:
            spaced = ''
            for l in range(len(line[0])):
                if (l+1)%3==0:
                    spaced+=str(line[0][l])+' '
                else:
                    spaced+=str(line[0][l])
            output_data+=str(spaced.strip('\n'))+','+str(line[1])+'\n'
            #output_data+=str(line[0]).strip('\n')+','+str(line[1]).strip('\n')+'\n'

        response = make_response(output_data)

        response.headers["Content-Disposition"] = "attachment; filename=High_Quality_Sequences_Round9.csv"
        return response

    if request.method == "GET":
        return render_template("quality.html")
@app.route('/download_qual_pre_9', methods=["GET", "POST"])
def download_qual_pre_9():
    if request.method == 'POST':

        if shared_dict["trim_type"] == 'pos':
            high_qual_seqs_not_unique, origonal_seqs_not_unique = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["num_start"], shared_dict["num_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))
        if shared_dict["trim_type"] == 'seq':
            high_qual_seqs_not_unique, origonal_seqs_not_unique = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["seq_start"], shared_dict["seq_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))

        rounds = range(len(shared_dict["all_high_qual_seqs"]))

        input_data = high_qual_seqs_not_unique[9]
        output_data = ''
        output_data+='Sequence,Reads\n'

        raw_data_list = [(k, v) for k, v in collections.Counter(input_data).items()]
        top_sequences = sorted(raw_data_list, key = lambda x: x[1], reverse=True)[:10000]

        for line in top_sequences:
            spaced = ''
            for l in range(len(line[0])):
                if (l+1)%3==0:
                    spaced+=str(line[0][l])+' '
                else:
                    spaced+=str(line[0][l])
            output_data+=str(spaced.strip('\n'))+','+str(line[1])+'\n'
            #output_data+=str(line[0]).strip('\n')+','+str(line[1]).strip('\n')+'\n'

        response = make_response(output_data)

        response.headers["Content-Disposition"] = "attachment; filename=High_Quality_Sequences_Round10.csv"
        return response

    if request.method == "GET":
        return render_template("quality.html")


##----------------- CONSOLIDATE W LOOP------------------##
@app.route('/download_qual_0', methods=["GET", "POST"])
def download_qual0():
    if request.method == 'POST':

        if shared_dict["trim_type"] == 'pos':
            high_qual_seqs_not_unique, origonal_seqs_not_unique = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["num_start"], shared_dict["num_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))
        if shared_dict["trim_type"] == 'seq':
            high_qual_seqs_not_unique, origonal_seqs_not_unique = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["seq_start"], shared_dict["seq_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))

        rounds = range(len(shared_dict["all_high_qual_seqs"]))

        input_data = high_qual_seqs_not_unique[0]
        output_data = ''
        output_data+='Sequence,Reads\n'

        raw_data_list = [(k, v) for k, v in collections.Counter(input_data).items()]
        top_sequences = sorted(raw_data_list, key = lambda x: x[1], reverse=True)

        for line in top_sequences:
            spaced = ''
            for l in range(len(line[0])):
                if (l+1)%3==0:
                    spaced+=str(line[0][l])+' '
                else:
                    spaced+=str(line[0][l])
            #output_data+=str(spaced).strip('\n')+','+str(line[1]).strip('\n')+'\n'
            output_data+=str(line[0]).strip('\n')+','+str(line[1]).strip('\n')+'\n'

        response = make_response(output_data)

        response.headers["Content-Disposition"] = "attachment; filename=High_Quality_Preprocessed_Sequences_Round1.csv"
        return response

    if request.method == "GET":
        return render_template("quality.html")

@app.route('/download_qual_1', methods=["GET", "POST"])
def download_qual1():
    if request.method == 'POST':

        if shared_dict["trim_type"] == 'pos':
            high_qual_seqs_not_unique, origonal_seqs_not_unique = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["num_start"], shared_dict["num_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))
        if shared_dict["trim_type"] == 'seq':
            high_qual_seqs_not_unique, origonal_seqs_not_unique = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["seq_start"], shared_dict["seq_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))

        rounds = range(len(shared_dict["all_high_qual_seqs"]))

        input_data = high_qual_seqs_not_unique[1]
        output_data = ''
        output_data+='Sequence,Reads\n'

        raw_data_list = [(k, v) for k, v in collections.Counter(input_data).items()]
        top_sequences = sorted(raw_data_list, key = lambda x: x[1], reverse=True)[:10000]

        for line in top_sequences:
            spaced = ''
            for l in range(len(line[0])):
                if (l+1)%3==0:
                    spaced+=str(line[0][l])+' '
                else:
                    spaced+=str(line[0][l])
            #output_data+=str(spaced.strip('\n'))+','+str(line[1])+'\n'
            output_data+=str(line[0]).strip('\n')+','+str(line[1]).strip('\n')+'\n'

        response = make_response(output_data)

        response.headers["Content-Disposition"] = "attachment; filename=High_Quality_Preprocessed_Sequences_Round2.csv"
        return response

    if request.method == "GET":
        return render_template("quality.html")
@app.route('/download_qual_2', methods=["GET", "POST"])
def download_qual2():
    if request.method == 'POST':

        if shared_dict["trim_type"] == 'pos':
            high_qual_seqs_not_unique, origonal_seqs_not_unique = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["num_start"], shared_dict["num_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))
        if shared_dict["trim_type"] == 'seq':
            high_qual_seqs_not_unique, origonal_seqs_not_unique = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["seq_start"], shared_dict["seq_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))

        rounds = range(len(shared_dict["all_high_qual_seqs"]))

        input_data = high_qual_seqs_not_unique[2]
        output_data = ''
        output_data+='Sequence,Reads\n'

        raw_data_list = [(k, v) for k, v in collections.Counter(input_data).items()]
        top_sequences = sorted(raw_data_list, key = lambda x: x[1], reverse=True)[:10000]

        for line in top_sequences:
            spaced = ''
            for l in range(len(line[0])):
                if (l+1)%3==0:
                    spaced+=str(line[0][l])+' '
                else:
                    spaced+=str(line[0][l])
            #output_data+=str(spaced.strip('\n'))+','+str(line[1])+'\n'
            output_data+=str(line[0]).strip('\n')+','+str(line[1]).strip('\n')+'\n'

        response = make_response(output_data)

        response.headers["Content-Disposition"] = "attachment; filename=High_Quality_Preprocessed_Sequences_Round3.csv"
        return response

    if request.method == "GET":
        return render_template("quality.html")
@app.route('/download_qual_3', methods=["GET", "POST"])
def download_qual3():
    if request.method == 'POST':

        if shared_dict["trim_type"] == 'pos':
            high_qual_seqs_not_unique, origonal_seqs_not_unique = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["num_start"], shared_dict["num_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))
        if shared_dict["trim_type"] == 'seq':
            high_qual_seqs_not_unique, origonal_seqs_not_unique = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["seq_start"], shared_dict["seq_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))

        rounds = range(len(shared_dict["all_high_qual_seqs"]))

        input_data = high_qual_seqs_not_unique[3]
        output_data = ''
        output_data+='Sequence,Reads\n'

        raw_data_list = [(k, v) for k, v in collections.Counter(input_data).items()]
        top_sequences = sorted(raw_data_list, key = lambda x: x[1], reverse=True)[:10000]

        for line in top_sequences:
            spaced = ''
            for l in range(len(line[0])):
                if (l+1)%3==0:
                    spaced+=str(line[0][l])+' '
                else:
                    spaced+=str(line[0][l])
            #output_data+=str(spaced.strip('\n'))+','+str(line[1])+'\n'
            output_data+=str(line[0]).strip('\n')+','+str(line[1]).strip('\n')+'\n'

        response = make_response(output_data)

        response.headers["Content-Disposition"] = "attachment; filename=High_Quality_Preprocessed_Sequences_Round4.csv"
        return response

    if request.method == "GET":
        return render_template("quality.html")
@app.route('/download_qual_4', methods=["GET", "POST"])
def download_qual4():
    if request.method == 'POST':

        if shared_dict["trim_type"] == 'pos':
            high_qual_seqs_not_unique, origonal_seqs_not_unique = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["num_start"], shared_dict["num_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))
        if shared_dict["trim_type"] == 'seq':
            high_qual_seqs_not_unique, origonal_seqs_not_unique = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["seq_start"], shared_dict["seq_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))

        rounds = range(len(shared_dict["all_high_qual_seqs"]))

        input_data = high_qual_seqs_not_unique[4]
        output_data = ''
        output_data+='Sequence,Reads\n'

        raw_data_list = [(k, v) for k, v in collections.Counter(input_data).items()]
        top_sequences = sorted(raw_data_list, key = lambda x: x[1], reverse=True)[:10000]

        for line in top_sequences:
            spaced = ''
            for l in range(len(line[0])):
                if (l+1)%3==0:
                    spaced+=str(line[0][l])+' '
                else:
                    spaced+=str(line[0][l])
            #output_data+=str(spaced.strip('\n'))+','+str(line[1])+'\n'
            output_data+=str(line[0]).strip('\n')+','+str(line[1]).strip('\n')+'\n'

        response = make_response(output_data)

        response.headers["Content-Disposition"] = "attachment; filename=High_Quality_Preprocessed_Sequences_Round5.csv"
        return response

    if request.method == "GET":
        return render_template("quality.html")
@app.route('/download_qual_5', methods=["GET", "POST"])
def download_qual5():
    if request.method == 'POST':

        rounds = range(len(shared_dict["all_high_qual_seqs"]))

        if shared_dict["trim_type"] == 'pos':
            high_qual_seqs_not_unique, origonal_seqs_not_unique = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["num_start"], shared_dict["num_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))
        if shared_dict["trim_type"] == 'seq':
            high_qual_seqs_not_unique, origonal_seqs_not_unique = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["seq_start"], shared_dict["seq_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))

        rounds = range(len(shared_dict["all_high_qual_seqs"]))

        input_data = high_qual_seqs_not_unique[5]
        output_data = ''
        output_data+='Sequence,Reads\n'

        raw_data_list = [(k, v) for k, v in collections.Counter(input_data).items()]
        top_sequences = sorted(raw_data_list, key = lambda x: x[1], reverse=True)[:10000]

        for line in top_sequences:
            spaced = ''
            for l in range(len(line[0])):
                if (l+1)%3==0:
                    spaced+=str(line[0][l])+' '
                else:
                    spaced+=str(line[0][l])
            output_data+=str(spaced.strip('\n'))+','+str(line[1])+'\n'
            #output_data+=str(line[0]).strip('\n')+','+str(line[1]).strip('\n')+'\n'

        response = make_response(output_data)

        response.headers["Content-Disposition"] = "attachment; filename=High_Quality_Preprocessed_Sequences_Round6.csv"
        return response

    if request.method == "GET":
        return render_template("quality.html")
@app.route('/download_qual_6', methods=["GET", "POST"])
def download_qual6():
    if request.method == 'POST':

        if shared_dict["trim_type"] == 'pos':
            high_qual_seqs_not_unique, origonal_seqs_not_unique = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["num_start"], shared_dict["num_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))
        if shared_dict["trim_type"] == 'seq':
            high_qual_seqs_not_unique, origonal_seqs_not_unique = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["seq_start"], shared_dict["seq_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))

        rounds = range(len(shared_dict["all_high_qual_seqs"]))

        input_data = high_qual_seqs_not_unique[6]
        output_data = ''
        output_data+='Sequence,Reads\n'

        raw_data_list = [(k, v) for k, v in collections.Counter(input_data).items()]
        top_sequences = sorted(raw_data_list, key = lambda x: x[1], reverse=True)[:10000]

        for line in top_sequences:
            spaced = ''
            for l in range(len(line[0])):
                if (l+1)%3==0:
                    spaced+=str(line[0][l])+' '
                else:
                    spaced+=str(line[0][l])
            #output_data+=str(spaced.strip('\n'))+','+str(line[1])+'\n'
            output_data+=str(line[0]).strip('\n')+','+str(line[1]).strip('\n')+'\n'

        response = make_response(output_data)

        response.headers["Content-Disposition"] = "attachment; filename=High_Quality_Preprocessed_Sequences_Round7.csv"
        return response

    if request.method == "GET":
        return render_template("quality.html")
@app.route('/download_qual_7', methods=["GET", "POST"])
def download_qual7():
    if request.method == 'POST':

        rounds = range(len(shared_dict["all_high_qual_seqs"]))

        if shared_dict["trim_type"] == 'pos':
            high_qual_seqs_not_unique, origonal_seqs_not_unique = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["num_start"], shared_dict["num_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))
        if shared_dict["trim_type"] == 'seq':
            high_qual_seqs_not_unique, origonal_seqs_not_unique = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["seq_start"], shared_dict["seq_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))

        rounds = range(len(shared_dict["all_high_qual_seqs"]))

        input_data = high_qual_seqs_not_unique[7]
        output_data = ''
        output_data+='Sequence,Reads\n'

        raw_data_list = [(k, v) for k, v in collections.Counter(input_data).items()]
        top_sequences = sorted(raw_data_list, key = lambda x: x[1], reverse=True)[:10000]

        for line in top_sequences:
            spaced = ''
            for l in range(len(line[0])):
                if (l+1)%3==0:
                    spaced+=str(line[0][l])+' '
                else:
                    spaced+=str(line[0][l])
            output_data+=str(spaced.strip('\n'))+','+str(line[1])+'\n'
            #output_data+=str(line[0]).strip('\n')+','+str(line[1]).strip('\n')+'\n'

        response = make_response(output_data)

        response.headers["Content-Disposition"] = "attachment; filename=High_Quality_Preprocessed_Sequences_Round8.csv"
        return response

    if request.method == "GET":
        return render_template("quality.html")
@app.route('/download_qual_8', methods=["GET", "POST"])
def download_qual8():
    if request.method == 'POST':

        if shared_dict["trim_type"] == 'pos':
            high_qual_seqs_not_unique, origonal_seqs_not_unique = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["num_start"], shared_dict["num_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))
        if shared_dict["trim_type"] == 'seq':
            high_qual_seqs_not_unique, origonal_seqs_not_unique = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["seq_start"], shared_dict["seq_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))

        rounds = range(len(shared_dict["all_high_qual_seqs"]))

        input_data = high_qual_seqs_not_unique[8]
        output_data = ''
        output_data+='Sequence,Reads\n'

        raw_data_list = [(k, v) for k, v in collections.Counter(input_data).items()]
        top_sequences = sorted(raw_data_list, key = lambda x: x[1], reverse=True)[:10000]

        for line in top_sequences:
            spaced = ''
            for l in range(len(line[0])):
                if (l+1)%3==0:
                    spaced+=str(line[0][l])+' '
                else:
                    spaced+=str(line[0][l])
            output_data+=str(spaced.strip('\n'))+','+str(line[1])+'\n'
            #output_data+=str(line[0]).strip('\n')+','+str(line[1]).strip('\n')+'\n'

        response = make_response(output_data)

        response.headers["Content-Disposition"] = "attachment; filename=High_Quality_Preprocessed_Sequences_Round9.csv"
        return response

    if request.method == "GET":
        return render_template("quality.html")
@app.route('/download_qual_9', methods=["GET", "POST"])
def download_qual9():
    if request.method == 'POST':

        if shared_dict["trim_type"] == 'pos':
            high_qual_seqs_not_unique, origonal_seqs_not_unique = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["num_start"], shared_dict["num_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))
        if shared_dict["trim_type"] == 'seq':
            high_qual_seqs_not_unique, origonal_seqs_not_unique = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["seq_start"], shared_dict["seq_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))

        rounds = range(len(shared_dict["all_high_qual_seqs"]))

        input_data = high_qual_seqs_not_unique[9]
        output_data = ''
        output_data+='Sequence,Reads\n'

        raw_data_list = [(k, v) for k, v in collections.Counter(input_data).items()]
        top_sequences = sorted(raw_data_list, key = lambda x: x[1], reverse=True)[:10000]

        for line in top_sequences:
            spaced = ''
            for l in range(len(line[0])):
                if (l+1)%3==0:
                    spaced+=str(line[0][l])+' '
                else:
                    spaced+=str(line[0][l])
            output_data+=str(spaced.strip('\n'))+','+str(line[1])+'\n'
            #output_data+=str(line[0]).strip('\n')+','+str(line[1]).strip('\n')+'\n'

        response = make_response(output_data)

        response.headers["Content-Disposition"] = "attachment; filename=High_Quality_Preprocessed_Sequences_Round10.csv"
        return response

    if request.method == "GET":
        return render_template("quality.html")

####------------------------------_########
@app.route('/download_clusters_0', methods=["GET", "POST"])
def download_clusters0():
    if request.method == 'POST':

        round_num = 0
        rounds = range(len(shared_dict["all_high_qual_seqs"]))
        if shared_dict["trim_type"] == 'pos':
            high_qual_seqs_not_unique, origonal_seqs_not_unique = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["num_start"], shared_dict["num_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))
        if shared_dict["trim_type"] == 'seq':
            high_qual_seqs_not_unique, origonal_seqs_not_unique = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["seq_start"], shared_dict["seq_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))

        input_data = list(set(shared_dict["all_high_qual_seqs"][round_num]))
        output_data = ''
        output_data+='Sequence,Cluster,Reads\n'
        used = []
        range_max = min(len(shared_dict["all_labels"][round_num]), 1000)
        for line in range(range_max):
            if input_data[line] not in used:
                used.append(input_data[line])
                cluster = shared_dict["all_labels"][round_num][line]
                counts = origonal_seqs_not_unique[round_num].count(input_data[line])
                #output_data+=str(input_data[line].strip('\n'))+','+str(cluster)+','+str(shared_dict["origonal_seqs"][round_num].count(input_data[line+"\n"]))+'\n'
                output_data+=str(input_data[line].strip('\n'))+','+str(cluster)+','+str(counts)+'\n'
        response = make_response(output_data)
        response.headers["Content-Disposition"] = "attachment; filename=Sequence_Clusters.csv"
        return response

    if request.method == "GET":
        return render_template("cluster_seqs.html")

@app.route('/download_clusters_1', methods=["GET", "POST"])
def download_clusters1():
    if request.method == 'POST':

        round_num = 1
        rounds = range(len(shared_dict["all_high_qual_seqs"]))
        if shared_dict["trim_type"] == 'pos':
            high_qual_seqs_not_unique, origonal_seqs_not_unique = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["num_start"], shared_dict["num_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))
        if shared_dict["trim_type"] == 'seq':
            high_qual_seqs_not_unique, origonal_seqs_not_unique = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["seq_start"], shared_dict["seq_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))

        input_data = shared_dict["all_high_qual_seqs"][round_num]
        output_data = ''
        output_data+='Sequence,Cluster,Reads\n'
        used = []
        range_max = min(len(shared_dict["all_labels"][round_num]), 1000)
        for line in range(range_max):
            if input_data[line] not in used:
                used.append(input_data[line])
                cluster = shared_dict["all_labels"][round_num][line]
                counts = origonal_seqs_not_unique[round_num].count(input_data[line])
                #output_data+=str(input_data[line].strip('\n'))+','+str(cluster)+','+str(shared_dict["origonal_seqs"][round_num].count(input_data[line+"\n"]))+'\n'
                output_data+=str(input_data[line].strip('\n'))+','+str(cluster)+','+str(counts)+'\n'
        response = make_response(output_data)
        response.headers["Content-Disposition"] = "attachment; filename=Sequence_Clusters.csv"
        return response

    if request.method == "GET":
        return render_template("cluster_seqs.html")

@app.route('/download_clusters_2', methods=["GET", "POST"])
def download_clusters2():
    if request.method == 'POST':

        round_num = 2
        rounds = range(len(shared_dict["all_high_qual_seqs"]))
        if shared_dict["trim_type"] == 'pos':
            high_qual_seqs_not_unique, origonal_seqs_not_unique = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["num_start"], shared_dict["num_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))
        if shared_dict["trim_type"] == 'seq':
            high_qual_seqs_not_unique, origonal_seqs_not_unique = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["seq_start"], shared_dict["seq_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))

        input_data = shared_dict["all_high_qual_seqs"][round_num]
        output_data = ''
        output_data+='Sequence,Cluster,Reads\n'
        used = []
        range_max = min(len(shared_dict["all_labels"][round_num]), 1000)
        for line in range(range_max):
            if input_data[line] not in used:
                used.append(input_data[line])
                cluster = shared_dict["all_labels"][round_num][line]
                counts = origonal_seqs_not_unique[round_num].count(input_data[line])
                #output_data+=str(input_data[line].strip('\n'))+','+str(cluster)+','+str(shared_dict["origonal_seqs"][round_num].count(input_data[line+"\n"]))+'\n'
                output_data+=str(input_data[line].strip('\n'))+','+str(cluster)+','+str(counts)+'\n'
        response = make_response(output_data)
        response.headers["Content-Disposition"] = "attachment; filename=Sequence_Clusters.csv"
        return response
    if request.method == "GET":
        return render_template("cluster_seqs.html")

@app.route('/download_clusters_3', methods=["GET", "POST"])
def download_clusters3():
    if request.method == 'POST':

        round_num = 3
        rounds = range(len(shared_dict["all_high_qual_seqs"]))
        if shared_dict["trim_type"] == 'pos':
            high_qual_seqs_not_unique, origonal_seqs_not_unique = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["num_start"], shared_dict["num_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))
        if shared_dict["trim_type"] == 'seq':
            high_qual_seqs_not_unique, origonal_seqs_not_unique = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["seq_start"], shared_dict["seq_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))

        input_data = shared_dict["all_high_qual_seqs"][round_num]
        output_data = ''
        output_data+='Sequence,Cluster,Reads\n'
        used = []
        range_max = min(len(shared_dict["all_labels"][round_num]), 1000)
        for line in range(range_max):
            if input_data[line] not in used:
                used.append(input_data[line])
                cluster = shared_dict["all_labels"][round_num][line]
                counts = origonal_seqs_not_unique[round_num].count(input_data[line])
                #output_data+=str(input_data[line].strip('\n'))+','+str(cluster)+','+str(shared_dict["origonal_seqs"][round_num].count(input_data[line+"\n"]))+'\n'
                output_data+=str(input_data[line].strip('\n'))+','+str(cluster)+','+str(counts)+'\n'
        response = make_response(output_data)
        response.headers["Content-Disposition"] = "attachment; filename=Sequence_Clusters.csv"
        return response

    if request.method == "GET":
        return render_template("cluster_seqs.html")

@app.route('/download_clusters_4', methods=["GET", "POST"])
def download_clusters4():
    if request.method == 'POST':

        round_num = 4
        rounds = range(len(shared_dict["all_high_qual_seqs"]))

        if shared_dict["trim_type"] == 'pos':
            high_qual_seqs_not_unique, origonal_seqs_not_unique = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["num_start"], shared_dict["num_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))
        if shared_dict["trim_type"] == 'seq':
            high_qual_seqs_not_unique, origonal_seqs_not_unique = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["seq_start"], shared_dict["seq_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))

        input_data = shared_dict["all_high_qual_seqs"][round_num]
        output_data = ''
        output_data+='Sequence,Cluster,Reads\n'
        used = []
        range_max = min(len(shared_dict["all_labels"][round_num]), 1000)
        for line in range(range_max):
            if input_data[line] not in used:
                used.append(input_data[line])
                cluster = shared_dict["all_labels"][round_num][line]
                counts = origonal_seqs_not_unique[round_num].count(input_data[line])
                #output_data+=str(input_data[line].strip('\n'))+','+str(cluster)+','+str(shared_dict["origonal_seqs"][round_num].count(input_data[line+"\n"]))+'\n'
                output_data+=str(input_data[line].strip('\n'))+','+str(cluster)+','+str(counts)+'\n'
        response = make_response(output_data)
        response.headers["Content-Disposition"] = "attachment; filename=Sequence_Clusters.csv"
        return response

    if request.method == "GET":
        return render_template("cluster_seqs.html")

@app.route('/download_clusters_5', methods=["GET", "POST"])
def download_clusters5():
    if request.method == 'POST':

        round_num = 4
        rounds = range(len(shared_dict["all_high_qual_seqs"]))

        if shared_dict["trim_type"] == 'pos':
            high_qual_seqs_not_unique, origonal_seqs_not_unique = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["num_start"], shared_dict["num_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))
        if shared_dict["trim_type"] == 'seq':
            high_qual_seqs_not_unique, origonal_seqs_not_unique = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["seq_start"], shared_dict["seq_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))

        input_data = shared_dict["all_high_qual_seqs"][round_num]
        output_data = ''
        output_data+='Sequence,Cluster,Reads\n'
        used = []
        range_max = min(len(shared_dict["all_labels"][round_num]), 1000)
        for line in range(range_max):
            if input_data[line] not in used:
                used.append(input_data[line])
                cluster = shared_dict["all_labels"][round_num][line]
                counts = origonal_seqs_not_unique[round_num].count(input_data[line])
                #output_data+=str(input_data[line].strip('\n'))+','+str(cluster)+','+str(shared_dict["origonal_seqs"][round_num].count(input_data[line+"\n"]))+'\n'
                output_data+=str(input_data[line].strip('\n'))+','+str(cluster)+','+str(counts)+'\n'
        response = make_response(output_data)
        response.headers["Content-Disposition"] = "attachment; filename=Sequence_Clusters.csv"
        return response

    if request.method == "GET":
        return render_template("cluster_seqs.html")

@app.route('/download_clusters_6', methods=["GET", "POST"])
def download_clusters6():
    if request.method == 'POST':

        round_num = 6
        rounds = range(len(shared_dict["all_high_qual_seqs"]))
        if shared_dict["trim_type"] == 'pos':
            high_qual_seqs_not_unique, origonal_seqs_not_unique = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["num_start"], shared_dict["num_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))
        if shared_dict["trim_type"] == 'seq':
            high_qual_seqs_not_unique, origonal_seqs_not_unique = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["seq_start"], shared_dict["seq_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))


        input_data = shared_dict["all_high_qual_seqs"][round_num]
        output_data = ''
        output_data+='Sequence,Cluster,Reads\n'
        used = []
        range_max = min(len(shared_dict["all_labels"][round_num]), 1000)
        for line in range(range_max):
            if input_data[line] not in used:
                used.append(input_data[line])
                cluster = shared_dict["all_labels"][round_num][line]
                counts = origonal_seqs_not_unique[round_num].count(input_data[line])
                #output_data+=str(input_data[line].strip('\n'))+','+str(cluster)+','+str(shared_dict["origonal_seqs"][round_num].count(input_data[line+"\n"]))+'\n'
                output_data+=str(input_data[line].strip('\n'))+','+str(cluster)+','+str(counts)+'\n'
        response = make_response(output_data)
        response.headers["Content-Disposition"] = "attachment; filename=Sequence_Clusters.csv"
        return response

    if request.method == "GET":
        return render_template("cluster_seqs.html")

@app.route('/download_clusters_7', methods=["GET", "POST"])
def download_clusters7():
    if request.method == 'POST':

        round_num = 7
        rounds = range(len(shared_dict["all_high_qual_seqs"]))
        if shared_dict["trim_type"] == 'pos':
            high_qual_seqs_not_unique, origonal_seqs_not_unique = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["num_start"], shared_dict["num_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))
        if shared_dict["trim_type"] == 'seq':
            high_qual_seqs_not_unique, origonal_seqs_not_unique = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["seq_start"], shared_dict["seq_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))

        input_data = shared_dict["all_high_qual_seqs"][round_num]
        output_data = ''
        output_data+='Sequence,Cluster,Reads\n'
        used = []
        range_max = min(len(shared_dict["all_labels"][round_num]), 1000)
        for line in range(range_max):
            if input_data[line] not in used:
                used.append(input_data[line])
                cluster = shared_dict["all_labels"][round_num][line]
                counts = origonal_seqs_not_unique[round_num].count(input_data[line])
                #output_data+=str(input_data[line].strip('\n'))+','+str(cluster)+','+str(shared_dict["origonal_seqs"][round_num].count(input_data[line+"\n"]))+'\n'
                output_data+=str(input_data[line].strip('\n'))+','+str(cluster)+','+str(counts)+'\n'
        response = make_response(output_data)
        response.headers["Content-Disposition"] = "attachment; filename=Sequence_Clusters.csv"
        return response

    if request.method == "GET":
        return render_template("cluster_seqs.html")

@app.route('/download_clusters_8', methods=["GET", "POST"])
def download_clusters8():
    if request.method == 'POST':

        round_num = 8
        rounds = range(len(shared_dict["all_high_qual_seqs"]))
        if shared_dict["trim_type"] == 'pos':
            high_qual_seqs_not_unique, origonal_seqs_not_unique = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["num_start"], shared_dict["num_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))
        if shared_dict["trim_type"] == 'seq':
            high_qual_seqs_not_unique, origonal_seqs_not_unique = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["seq_start"], shared_dict["seq_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))

        input_data = shared_dict["all_high_qual_seqs"][round_num]
        output_data = ''
        output_data+='Sequence,Cluster,Reads\n'
        used = []
        range_max = min(len(shared_dict["all_labels"][round_num]), 1000)
        for line in range(range_max):
            if input_data[line] not in used:
                used.append(input_data[line])
                cluster = shared_dict["all_labels"][round_num][line]
                counts = origonal_seqs_not_unique[round_num].count(input_data[line])
                #output_data+=str(input_data[line].strip('\n'))+','+str(cluster)+','+str(shared_dict["origonal_seqs"][round_num].count(input_data[line+"\n"]))+'\n'
                output_data+=str(input_data[line].strip('\n'))+','+str(cluster)+','+str(counts)+'\n'
        response = make_response(output_data)
        response.headers["Content-Disposition"] = "attachment; filename=Sequence_Clusters.csv"
        return response

    if request.method == "GET":
        return render_template("cluster_seqs.html")

@app.route('/download_clusters_9', methods=["GET", "POST"])
def download_clusters9():
    if request.method == 'POST':

        round_num = 9
        rounds = range(len(shared_dict["all_high_qual_seqs"]))
        if shared_dict["trim_type"] == 'pos':
            high_qual_seqs_not_unique, origonal_seqs_not_unique = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["num_start"], shared_dict["num_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))
        if shared_dict["trim_type"] == 'seq':
            high_qual_seqs_not_unique, origonal_seqs_not_unique = get_high_qual_seqs(shared_dict["trim_type"], shared_dict["qual_yes"],shared_dict["uploaded_files"], shared_dict["seq_start"], shared_dict["seq_end"], shared_dict["rc"], shared_dict["quality_cutoff"], range(shared_dict["num_rounds"]))


        input_data = shared_dict["all_high_qual_seqs"][round_num]
        output_data = ''
        output_data+='Sequence,Cluster,Reads\n'
        used = []
        range_max = min(len(shared_dict["all_labels"][round_num]), 1000)
        for line in range(range_max):
            if input_data[line] not in used:
                used.append(input_data[line])
                cluster = shared_dict["all_labels"][round_num][line]
                counts = origonal_seqs_not_unique[round_num].count(input_data[line])
                #output_data+=str(input_data[line].strip('\n'))+','+str(cluster)+','+str(shared_dict["origonal_seqs"][round_num].count(input_data[line+"\n"]))+'\n'
                output_data+=str(input_data[line].strip('\n'))+','+str(cluster)+','+str(counts)+'\n'
        response = make_response(output_data)
        response.headers["Content-Disposition"] = "attachment; filename=Sequence_Clusters.csv"
        return response

    if request.method == "GET":
        return render_template("cluster_seqs.html")

#########--------------__##############
@app.route('/download_region_0', methods=["GET", "POST"])
def download_region0():
    if request.method == 'POST':
        round_num = 0
        input_data = shared_dict["all_regions_all_counts_rounds_seqs"][round_num][0]
        output_data = ''
        output_data+='Sequence,Reads\n'
        for line in range(min(1000, len(shared_dict["all_regions_all_counts_rounds_seqs"][round_num][0]))):
            output_data+=str(input_data[line])+','+str(shared_dict["origonal_seqs"][round_num].count(input_data[line]))+'\n'
            #output_data+=str(input_data[line].strip('\n'))+','+str(shared_dict["origonal_seqs"][round_num].count(input_data[line+"\n"]))+'\n'
        response = make_response(output_data)
        response.headers["Content-Disposition"] = "attachment; filename=Motif_Tracker_Round1.csv"
        return response
    if request.method == "GET":
        return render_template("region_emergence.html")
@app.route('/download_region_1', methods=["GET", "POST"])
def download_region1():
    if request.method == 'POST':
        round_num = 1
        input_data = shared_dict["all_regions_all_counts_rounds_seqs"][round_num][0]
        output_data = ''
        output_data+='Sequence,Reads\n'
        for line in range(min(1000, len(shared_dict["all_regions_all_counts_rounds_seqs"][round_num][0]))):
            output_data+=str(input_data[line])+','+str(shared_dict["origonal_seqs"][round_num].count(input_data[line]))+'\n'
            #output_data+=str(input_data[line].strip('\n'))+','+str(shared_dict["origonal_seqs"][round_num].count(input_data[line+"\n"]))+'\n'
        response = make_response(output_data)
        response.headers["Content-Disposition"] = "attachment; filename=Motif_Tracker_Round2.csv"
        return response
    if request.method == "GET":
        return render_template("region_emergence.html")
@app.route('/download_region_2', methods=["GET", "POST"])
def download_region2():
    if request.method == 'POST':
        round_num = 2
        input_data = shared_dict["all_regions_all_counts_rounds_seqs"][round_num][0]
        output_data = ''
        output_data+='Sequence,Reads\n'
        for line in range(min(1000, len(shared_dict["all_regions_all_counts_rounds_seqs"][round_num][0]))):
            output_data+=str(input_data[line])+','+str(shared_dict["origonal_seqs"][round_num].count(input_data[line]))+'\n'
            #output_data+=str(input_data[line].strip('\n'))+','+str(shared_dict["origonal_seqs"][round_num].count(input_data[line+"\n"]))+'\n'
        response = make_response(output_data)
        response.headers["Content-Disposition"] = "attachment; filename=Motif_Tracker_Round3.csv"
        return response
    if request.method == "GET":
        return render_template("region_emergence.html")
@app.route('/download_region_3', methods=["GET", "POST"])
def download_region3():
    if request.method == 'POST':
        round_num = 3
        input_data = shared_dict["all_regions_all_counts_rounds_seqs"][round_num][0]
        output_data = ''
        output_data+='Sequence,Reads\n'
        for line in range(min(1000, len(shared_dict["all_regions_all_counts_rounds_seqs"][round_num][0]))):
            output_data+=str(input_data[line])+','+str(shared_dict["origonal_seqs"][round_num].count(input_data[line]))+'\n'
            #output_data+=str(input_data[line].strip('\n'))+','+str(shared_dict["origonal_seqs"][round_num].count(input_data[line+"\n"]))+'\n'
        response = make_response(output_data)
        response.headers["Content-Disposition"] = "attachment; filename=Motif_Tracker_Round4.csv"
        return response
    if request.method == "GET":
        return render_template("region_emergence.html")
@app.route('/download_region_4', methods=["GET", "POST"])
def download_region4():
    if request.method == 'POST':
        round_num = 4
        input_data = shared_dict["all_regions_all_counts_rounds_seqs"][round_num][0]
        output_data = ''
        output_data+='Sequence,Reads\n'
        for line in range(min(1000, len(shared_dict["all_regions_all_counts_rounds_seqs"][round_num][0]))):
            output_data+=str(input_data[line])+','+str(shared_dict["origonal_seqs"][round_num].count(input_data[line]))+'\n'
            #output_data+=str(input_data[line].strip('\n'))+','+str(shared_dict["origonal_seqs"][round_num].count(input_data[line+"\n"]))+'\n'
        response = make_response(output_data)
        response.headers["Content-Disposition"] = "attachment; filename=Motif_Tracker_Round5.csv"
        return response
    if request.method == "GET":
        return render_template("region_emergence.html")
@app.route('/download_region_5', methods=["GET", "POST"])
def download_region5():
    if request.method == 'POST':
        round_num = 5
        input_data = shared_dict["all_regions_all_counts_rounds_seqs"][round_num][0]
        output_data = ''
        output_data+='Sequence,Reads\n'
        for line in range(min(1000, len(shared_dict["all_regions_all_counts_rounds_seqs"][round_num][0]))):
            output_data+=str(input_data[line])+','+str(shared_dict["origonal_seqs"][round_num].count(input_data[line]))+'\n'
            #output_data+=str(input_data[line].strip('\n'))+','+str(shared_dict["origonal_seqs"][round_num].count(input_data[line+"\n"]))+'\n'
        response = make_response(output_data)
        response.headers["Content-Disposition"] = "attachment; filename=Motif_Tracker_Round6.csv"
        return response
    if request.method == "GET":
        return render_template("region_emergence.html")
@app.route('/download_region_6', methods=["GET", "POST"])
def download_region6():
    if request.method == 'POST':
        round_num = 6
        input_data = shared_dict["all_regions_all_counts_rounds_seqs"][round_num][0]
        output_data = ''
        output_data+='Sequence,Reads\n'
        for line in range(min(1000, len(shared_dict["all_regions_all_counts_rounds_seqs"][round_num][0]))):
            output_data+=str(input_data[line])+','+str(shared_dict["origonal_seqs"][round_num].count(input_data[line]))+'\n'
            #output_data+=str(input_data[line].strip('\n'))+','+str(shared_dict["origonal_seqs"][round_num].count(input_data[line+"\n"]))+'\n'
        response = make_response(output_data)
        response.headers["Content-Disposition"] = "attachment; filename=Motif_Tracker_Round7.csv"
        return response
    if request.method == "GET":
        return render_template("region_emergence.html")
@app.route('/download_region_7', methods=["GET", "POST"])
def download_region7():
    if request.method == 'POST':
        round_num = 7
        input_data = shared_dict["all_regions_all_counts_rounds_seqs"][round_num][0]
        output_data = ''
        output_data+='Sequence,Reads\n'
        for line in range(min(1000, len(shared_dict["all_regions_all_counts_rounds_seqs"][round_num][0]))):
            output_data+=str(input_data[line])+','+str(shared_dict["origonal_seqs"][round_num].count(input_data[line]))+'\n'
            #output_data+=str(input_data[line].strip('\n'))+','+str(shared_dict["origonal_seqs"][round_num].count(input_data[line+"\n"]))+'\n'
        response = make_response(output_data)
        response.headers["Content-Disposition"] = "attachment; filename=Motif_Tracker_Round8.csv"
        return response
    if request.method == "GET":
        return render_template("region_emergence.html")
@app.route('/download_region_8', methods=["GET", "POST"])
def download_region8():
    if request.method == 'POST':
        round_num = 8
        input_data = shared_dict["all_regions_all_counts_rounds_seqs"][round_num][0]
        output_data = ''
        output_data+='Sequence,Reads\n'
        for line in range(min(1000, len(shared_dict["all_regions_all_counts_rounds_seqs"][round_num][0]))):
            output_data+=str(input_data[line])+','+str(shared_dict["origonal_seqs"][round_num].count(input_data[line]))+'\n'
            #output_data+=str(input_data[line].strip('\n'))+','+str(shared_dict["origonal_seqs"][round_num].count(input_data[line+"\n"]))+'\n'
        response = make_response(output_data)
        response.headers["Content-Disposition"] = "attachment; filename=Motif_Tracker_Round9.csv"
        return response
    if request.method == "GET":
        return render_template("region_emergence.html")
@app.route('/download_region_9', methods=["GET", "POST"])
def download_region9():
    if request.method == 'POST':
        round_num = 9
        input_data = shared_dict["all_regions_all_counts_rounds_seqs"][round_num][0]
        output_data = ''
        output_data+='Sequence,Reads\n'
        for line in range(min(1000, len(shared_dict["all_regions_all_counts_rounds_seqs"][round_num][0]))):
            output_data+=str(input_data[line])+','+str(shared_dict["origonal_seqs"][round_num].count(input_data[line]))+'\n'
            #output_data+=str(input_data[line].strip('\n'))+','+str(shared_dict["origonal_seqs"][round_num].count(input_data[line+"\n"]))+'\n'
        response = make_response(output_data)
        response.headers["Content-Disposition"] = "attachment; filename=Motif_Tracker_Round10.csv"
        return response
    if request.method == "GET":
        return render_template("region_emergence.html")

if __name__ == "__main__":
    app.run()
