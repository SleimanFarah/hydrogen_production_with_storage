import pandas as pd

# Create a DataFrame for the tasks
data = {
    'Task': ['T1: Implement the system', 'T2: Monitor the system\'s performance', 'T3: Utilize feedback loops'],
    'Start Date': ['2024-06-03', '2024-07-01', '2024-09-01'],
    'Duration (weeks)': [4, 8, 6]
}

df = pd.DataFrame(data)

# Convert Start Date to datetime
df['Start Date'] = pd.to_datetime(df['Start Date'])

# Calculate End Date
df['End Date'] = df['Start Date'] + pd.to_timedelta(df['Duration (weeks)'], unit='W')

# Calculate Start (Days) from the project start date (2024-06-03)
project_start_date = pd.to_datetime('2024-06-03')
df['Start (Days)'] = (df['Start Date'] - project_start_date).dt.days

# Calculate Duration (Days)
df['Duration (Days)'] = df['Duration (weeks)'] * 7

# Select only the necessary columns for the Gantt chart
gantt_chart_df = df[['Task', 'Start Date', 'Duration (weeks)', 'End Date', 'Start (Days)', 'Duration (Days)']]

# Save DataFrame to Excel
file_path = 'WP3_Gantt_Chart_Colored_Bars.xlsx'
gantt_chart_df.to_excel(file_path, index=False)

file_path
