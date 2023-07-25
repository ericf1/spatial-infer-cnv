import pandas as pd
from find_normal import find_normal
from find_normal_unlabeled import find_normal_unlabeled
from find_normal_more import find_normal_more
from unannotated_find_normal import unannotated_find_normal
import traceback

df = pd.read_csv("specific_df.csv")

# get the AnnotationFile, LibraryName, and Spaceranger
df = df[["LibraryName", "AnnotationFile",  "Spaceranger"]]

# only do these libraries
# these_libraries = ['HT339B2-H1Fp1U1Bp1', 'HT397B1-S1H3A1U1Bp1', 'HT268B1-Th1H3Fc2U12Z1Bs1', 'HT268B1-Th1H3Fc2U22Z1Bs1', 'HT268B1-Th1H3Fc2U2Z1Bs1', 'HT268B1-Th1H3Fc2U32Z1Bs1', 'HT413C1-Th1K2A4U14Bp1', 'HT413C1-Th1K2A4U2Bp1',
#                    'HT448C1-Th1K1Fp1U13Bp1', 'HT448C1-Th1K1Fp1U1Bp1', 'HT397B1-S1H3A1U21Bp1', 'HT565B1-S1H2Fp2U1Bp1', 'HT472C1-S1H4Fp1U1Bp1', 'HT472C1-Th1K1Fp1U1Bp1', 'HT413C1-Th1K2A4U1Bp1', 'HT553P1-S1H2Fp1U1Bp1']

#library_list = ['HT480B1-S1H1Fs1U1Bp1']

#library_list = ['HT297B1-S1H1Fc2U1Z1Bs1', 'HT323B1-S1H1Fc2U1Z1Bs1', 'HT397B1-S1H3Fs1U1Bp1']
#library_list = ['HT206B1-U1_ST_Bn1']

# string = "HT206B1-S1H8_U1.csv HT206B1-S1H8_U3.csv HT206B1-S1H8_U5.csv HT308B1-S1H1.csv HT308B1-S1H5.csv HT339B1-S1H3_U1.csv HT339B2-S1H1.csv HT397B1-S1H3.csv HT480B1-S1H1_U2.csv HT206B1-S1H8_U2.csv HT206B1-S1H8_U4.csv HT297B1-S1H1.csv HT308B1-S1H4.csv HT323B1-S1H1.csv HT339B1-S1H3_U2.csv HT397B1-S1H2.csv HT480B1-S1H1_U1.csv HT565B1-S1H2.csv"

# Split the string into a list
# file_list = string.split()

# Remove the ".csv" extension from each element
# library_list = [filename.replace('.csv', '') for filename in file_list]

# print(library_list)

df = df[df["LibraryName"].isin(library_list)]

failed = []
for r in df.itertuples():
    try:
        print(r.LibraryName, "start")
        find_normal_unlabeled(
            r.Spaceranger,
            r.AnnotationFile,
            r.LibraryName,
            directory="/diskmnt/Projects/Users/efang/normalized_annotations_unlabeled/"
        )
        print(r.LibraryName, "done")
    except:
        traceback.print_exc()
        print(r.LibraryName, "failed")
        failed.append(r.LibraryName)
print(failed)
