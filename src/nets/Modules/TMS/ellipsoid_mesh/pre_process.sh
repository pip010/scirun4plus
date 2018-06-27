# MakeSoloMaterial.py

# cmmd = r'"%s" pad -min -4 -4 -4 -max M+4 M+4 M+4 -b pad -v 0 -i "%s" -o "%s"' % (unu_path, innrrd, outnrrd)
# unorient_cmd = r'"%s" -input "%s" -output "%s" -transform "%s"' % (os.path.join(binary_path,"UnorientNrrdAndGetTransform"), innrrd, unorient_fname, transform_fname)
# cmmd = r'"%s" make -i %s.lut.raw -t float -s %d -e ascii -o %s.lut.nrrd' % (unu_path, name, maxval+1, name)
# cmmd = r'"%s" lut -m %s.lut.nrrd -min 0 -max %d -t float -i "%s" -o "%s"' % (unu_path, name, maxval, innrrd, oname)
# cmmd = r'"%s" %s "%s" "%s"' % (os.path.join(binary_path,"morphsmooth"), r, innrd, oname)
# cmmd = "%s %s" % (os.path.join(binary_path,"ComputeTightenedLabels"), tightened_files)

#unu pad -min -4 -4 -4 -max M+4 M+4 M+4 -b pad -v 0 -i $1 | UnorientNrrdAndGetTransform -input - -output segmentation_unorient.nrrd -transform segmentation_transform.tf
UnorientNrrdAndGetTransform -input $1 -output segmentation_unorient.nrrd -transform segmentation_transform.tf

echo 1 0 0 0 0 | unu make -t float -s 5 -e ascii | unu lut -i segmentation_unorient.nrrd -m - -min 1 -max 6 -t float -o tmp_c0.solo.nrrd
echo 0 1 0 0 0 | unu make -t float -s 5 -e ascii | unu lut -i segmentation_unorient.nrrd -m - -min 1 -max 6 -t float -o tmp_c1.solo.nrrd
echo 0 0 1 0 0 | unu make -t float -s 5 -e ascii | unu lut -i segmentation_unorient.nrrd -m - -min 1 -max 6 -t float -o tmp_c2.solo.nrrd 
echo 0 0 0 1 0 | unu make -t float -s 5 -e ascii | unu lut -i segmentation_unorient.nrrd -m - -min 1 -max 6 -t float -o tmp_c3.solo.nrrd
echo 0 0 0 0 1 | unu make -t float -s 5 -e ascii | unu lut -i segmentation_unorient.nrrd -m - -min 1 -max 6 -t float -o tmp_c4.solo.nrrd


#ComputeTightenedLabels c1.tight.nrrd c2.tight.nrrd c3.tight.nrrd c4.tight.nrrd c5.tight.nrrd
morphsmooth $2 tmp_c0.solo.nrrd tmp_c0.tight.nrrd > tmp_morph0.txt &
morphsmooth $2 tmp_c1.solo.nrrd tmp_c1.tight.nrrd > tmp_morph1.txt &
morphsmooth $2 tmp_c2.solo.nrrd tmp_c2.tight.nrrd > tmp_morph2.txt &
morphsmooth $2 tmp_c3.solo.nrrd tmp_c3.tight.nrrd > tmp_morph3.txt &
morphsmooth $2 tmp_c4.solo.nrrd tmp_c4.tight.nrrd > tmp_morph4.txt &
wait

# Transform back

ConvertNrrdToField tmp_c0.tight.nrrd tmp_c0.fld
ConvertNrrdToField tmp_c1.tight.nrrd tmp_c1.fld
ConvertNrrdToField tmp_c2.tight.nrrd tmp_c2.fld
ConvertNrrdToField tmp_c3.tight.nrrd tmp_c3.fld
ConvertNrrdToField tmp_c4.tight.nrrd tmp_c4.fld

#ExtractIsosurface tmp_c0.fld tmp_c0_iso.ts
#ExtractIsosurface tmp_c1.fld tmp_c1_iso.ts
#ExtractIsosurface tmp_c2.fld tmp_c2_iso.ts
#ExtractIsosurface tmp_c3.fld tmp_c3_iso.ts
#ExtractIsosurface tmp_c4.fld tmp_c4_iso.ts
#ExtractIsosurface tmp_c5.fld tmp_c5_iso.ts

TransformFieldWithTransform -input tmp_c0.fld -transform segmentation_transform.tf -output tmp_c0_transformed.fld
TransformFieldWithTransform -input tmp_c1.fld -transform segmentation_transform.tf -output tmp_c1_transformed.fld
TransformFieldWithTransform -input tmp_c2.fld -transform segmentation_transform.tf -output tmp_c2_transformed.fld
TransformFieldWithTransform -input tmp_c3.fld -transform segmentation_transform.tf -output tmp_c3_transformed.fld
TransformFieldWithTransform -input tmp_c4.fld -transform segmentation_transform.tf -output tmp_c4_transformed.fld

ConvertFieldToNrrd tmp_c0_transformed.fld c0_transformed.nrrd
ConvertFieldToNrrd tmp_c1_transformed.fld c1_transformed.nrrd
ConvertFieldToNrrd tmp_c2_transformed.fld c2_transformed.nrrd
ConvertFieldToNrrd tmp_c3_transformed.fld c3_transformed.nrrd
ConvertFieldToNrrd tmp_c4_transformed.fld c4_transformed.nrrd

#compress NRRD
unu save -i c0_transformed.nrrd -f nrrd -e gzip -o c0_transformed.nrrd
unu save -i c1_transformed.nrrd -f nrrd -e gzip -o c1_transformed.nrrd
unu save -i c2_transformed.nrrd -f nrrd -e gzip -o c2_transformed.nrrd
unu save -i c3_transformed.nrrd -f nrrd -e gzip -o c3_transformed.nrrd
unu save -i c4_transformed.nrrd -f nrrd -e gzip -o c4_transformed.nrrd


# CLEAN tmps
rm tmp_*
