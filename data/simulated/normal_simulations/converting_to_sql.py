import sqlite3
import ast

with open('mph_mean_0.5_mph_std0.1_iph_mean0.5_iph_std_0.1.txt','r') as f:
    data = ast.literal_eval(f.read())

conn = sqlite3.connect('mph_mean_0.5_mph_std0.1_iph_mean0.5_iph_std_0.1.db')
c = conn.cursor()

# Create table
# c.execute('''
    # CREATE TABLE genes_phenotypes (
        # gene TEXT,
        # phenotype TEXT
    # )
# ''')

# conn.commit()


# Iterate over the data
# for item in data:
    # gene = item[0]
    # phenotypes = item[1]

    # # For each phenotype, insert a row into the table
    # for phenotype in phenotypes:
        # c.execute('''
            # INSERT INTO genes_phenotypes (gene, phenotype)
            # VALUES (?, ?)
        # ''', (gene, phenotype))

# Commit the changes and close the connection
# conn.commit()
# conn.close()


# Get all phenotypes for a specific gene
c.execute('''
    SELECT phenotype
    FROM genes_phenotypes
    WHERE gene = ?
''', ('16',))

# Fetch all the results
results = c.fetchall()

# Print the results
for row in results:
    print(row[0])

# Close the connection
conn.close()
