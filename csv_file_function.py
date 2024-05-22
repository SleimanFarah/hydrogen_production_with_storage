
import csv
with open('solar_capacity_factor.csv', 'r') as read_obj: # read csv file as a list of lists
  csv_reader = csv.reader(read_obj) # pass the file object to reader() to get the reader object
  list_of_rows = list(csv_reader) # Pass reader object to list() to get a list of lists

a = []

for i in range(1, len(list_of_rows)):
    a = a + [float(list_of_rows[i][1])]


print(list_of_rows)


