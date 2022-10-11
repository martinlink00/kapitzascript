{

rm params.txt
echo "3" > params.txt
math12 -script potential.wls
math12 -script incond.wls
rm params.txt
echo "2" > params.txt
math12 -script potential.wls
math12 -script incond.wls
rm params.txt
echo "4" > params.txt
math12 -script potential.wls
math12 -script incond.wls
math12 -script paramest.wls

# Try
} || {
# Executed when above fails
rm params.txt
echo "3" > params.txt
math -script potential.wls
math -script incond.wls
rm params.txt
echo "2" > params.txt
math -script potential.wls
math -script incond.wls
rm params.txt
echo "4" > params.txt
math -script potential.wls
math -script incond.wls
math -script paramest.wls

}
