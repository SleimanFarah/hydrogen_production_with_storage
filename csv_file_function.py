def csv_reader_function(read_obj):
    import csv # read csv file as a list of lists
    csv_reader = csv.reader(read_obj) # pass the file object to reader() to get the reader object
    list_of_rows = list(csv_reader) # Pass reader object to list() to get a list of lists

    a = []

    for i in range(1, len(list_of_rows)):
        a = a + [float(list_of_rows[i][1])]

    return a



