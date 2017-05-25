import pandas as pd
import requests
import json
import sys
import os
from bs4 import BeautifulSoup
from collections import OrderedDict

class Geno(object):
    """class to reads and interprets 23andMe raw genotype data."""
    def __init__(self, raw_23andMe):
        """input path to raw 23andMe genotype data. data will be read as a form of pandas DataFrame.

        >>> g = Geno("./23andMe_raw_data.txt")
        >>> g.my_geno"""
        # GET report as json from 23andMe API
        res = requests.get("https://api.23andme.com/3/report/").content
        self._report = json.loads(res)
        self.my_geno = pd.read_csv(raw_23andMe, sep="\t", header=14,
                                   names=["rsid", "chr", "pos", "genotype"],
                                   low_memory=False)

    def type_id(self, rsid, ref_population="JPT"):
        """search input rsid at SNPedia, report your possible trait related to this SNP.

        >>> g.type_id("rs28897696")"""
        mine = self.my_geno[self.my_geno.rsid == rsid]
        if mine.empty:
            return None

        # your genotype. if you are a male(have "Y" chromosome),
        # duplicate your X chromosome allele to make a two-lettered string
        if mine.chr.values[0] in ["X", "Y"]:
            if "Y" in self.my_geno.chr.unique():
                mm = mine.genotype.values[0]
                genotype = _sort_seq(mm + ";" + mm)
        else:
            mm, ff = mine.genotype.values[0]
            genotype = _sort_seq(mm + ";" + ff)

        # search this at snpedia
        snpedia_report, g_percent = _search_rsid(rsid, genotype, ref_population)
        snpedia_report = snpedia_report.strip()
        print rsid, "|", genotype, "|", snpedia_report, "|", "{}% of {}".format(g_percent, ref_population)

        return mine

    def type_trait(self, query, ref_population="JPT"):
        """search input query at SNPedia, report all related SNPs and your genotype.

        >>> g.type_trait("lung cancer")"""
        rsid_list = _search_anything(query)
        if not rsid_list[0].startswith("rs"):
            print "try {} instead".format(rsid_list)
            return

        for rsid in rsid_list:
            self.type_id(rsid, ref_population)

        detail_url = ("https://www.snpedia.com/index.php/" + query).replace(" ", "%20")
        print " * more details at", detail_url

    def report_wellness(self, ref_population="JPT"):
        """return wellness report provided by 23andMe API sevice.

        >>> g.report_wellness()"""
        try:
            height, width = map(int, os.popen('stty size', 'r').read().split())
        except:
            width = 100

        print "="*(width//2-9), "WELLNESS REPORT", "="*(width//2-9), "\n"

        for i in range(len(self._report["data"])):
            with_detail = pd.DataFrame(self._report["data"][i]).dropna(subset=["details"])
            if not "markers" in with_detail.index:
                continue
            for j in range(len(with_detail.ix["markers"].details)):
                #print with_detail.report_id#with_detail.ix["markers"].details[j][u'biological_explanation']
                rsid_report = with_detail.ix["markers"].details[j]["id"]
                variants_report = with_detail.ix["markers"].details[j]["variants"]
                for variant_dict in variants_report:
                    pos_report = variant_dict["end"]
                    pick = self.my_geno[self.my_geno.rsid == rsid_report]
                    if pick.empty:
                        continue

                    #print pick.pos.values[0], pos_report
                    #if pick.pos.values[0] == pos_report:
                    from_mom, from_dad = pick.genotype.values[0]
                    genotype = from_mom+";"+from_dad
                    snpedia_report, geno_percent = _search_rsid(rsid_report, genotype, ref_population)

                    if variant_dict["has_effect"] == True:
                        print "< {} >".format(with_detail.ix["markers"]["title"])
                        print "-" * width
                        print with_detail.ix["markers"].details[j]["biological_explanation"].strip()
                        print
                        print "  rsID:", rsid_report
                        print "  Yours:", genotype
                        if ("There is currently no text" not in snpedia_report) and (snpedia_report.strip() != ""):
                            print "  Detail:", snpedia_report.strip(), "|", "{}% of {}".format(geno_percent, ref_population)
                        print "\n\n"

def _search_rsid(rsid, genotype, ref_population="JPT"):
    """search input rsid at SNPedia, return and compare SNPedia result and yours"""
    genotype = _sort_seq(genotype)
    snpedia_query = "https://www.snpedia.com/index.php/" + rsid
    soup = BeautifulSoup(requests.get(snpedia_query).content, "lxml")
    table = soup.find_all("table", attrs={"class": "sortable"})
    stab_orien = soup.find_all("td")[3].text

    # complement my genotype if stabilized orientation is 'minus'
    if stab_orien == "minus":
        genotype = _comp_seq(genotype)

    if table != []:
        snpedia_var_dict = OrderedDict()
        for n, i in enumerate(table[0].find_all("td")):
            if n % 3 == 0: k = i.text.strip().strip("(").strip(")")
            if n % 3 == 1: m = i.text.strip()
            if n % 3 == 2:
                if len(k) == 2: continue
                snpedia_var_dict[k] = (m, i.text.strip())

        # save snpedia_report (result of input genotype)
        try:
            snpedia_report = " | ".join(snpedia_var_dict[genotype])
        except KeyError:
            snpedia_report = "No info about {}".format(genotype)

        # save geno_percent (which percent of people in JPT population have input genotype)
        try:
            geno_order = snpedia_var_dict.keys().index(genotype)
            #print snpedia_var_dict.keys(), genotype, rsid
            for td in soup.find_all("td"):
                if td.img is not None:
                    _chart = td.img["src"]
                    break

            _, vals, __, ___, pops, ____, _____, ______ = _chart.split("&")
            val1, val2, val3 = vals[6:].split("|")
            pops = pops[11:-2].split("|")

            try:
                jpt_idx = len(pops)-pops.index(ref_population)-1
                geno_ratio = [map(float, val1.split(","))[jpt_idx],
                              map(float, val2.split(","))[jpt_idx],
                              map(float, val3.split(","))[jpt_idx]]
                #print geno_ratio, geno_order
                geno_percent = geno_ratio[geno_order]
            except ValueError:
                geno_percent = "--"
        except ValueError:
            geno_percent = "--"
        except UnboundLocalError:
            geno_percent = "--"
    else:
        snpedia_report = ""
        geno_percent = "--"
    return snpedia_report, geno_percent

def _search_anything(query):
    """search input words at SNPedia, return and compare related genotypes and yours"""
    snpedia_query = "https://www.snpedia.com/index.php/" + query
    soup = BeautifulSoup(requests.get(snpedia_query).content, "lxml")
    links = soup.find_all("div", attrs={"class": "mw-content-ltr"})[0].find_all("a")
    rsid_list = [l.text for l in links if l.text.startswith("rs")]
    rsid_list = [txt for txt in rsid_list if not txt.endswith(")")]

    if rsid_list == []:
        rsid_list = [l.text for l in links if l.text]
        if 'search for this page title' in rsid_list:
            rsid_list = ["another phenotype"]
    return rsid_list

def _comp_seq(genotype):
    comp_dict = {"A":"T", "T":"A", "G":"C", "C":"G", ";":";"}
    return _sort_seq("".join([comp_dict[i] for i in genotype]))

def _sort_seq(genotype):
    order_dict = {"A":0, "C":1, "G":2, "T":3, ";":";"}
    m, _, f = [order_dict[i] for i in genotype]
    if m > f: return ";".join[f, m]
    else: return genotype
