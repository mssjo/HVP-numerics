# Splits the output of the hybrid plot mode into separate files for each integral

function separate(var) {
    target = "../plots/" file "_" var ".dat"
    for (i=2; i<NF; i++) { printf $i "\t" > target }
    print $NF > target
}

NR == 1 {
    separate("SJ1")
    separate("SJ2")
    separate("Hreg")
    separate("E1h2")
}
NR > 1 {
    separate($1)
}
