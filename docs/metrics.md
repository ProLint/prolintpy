# Contact-based metrics

Lipid-protein interactions are an association problem. During MD simulations, lipids (and proteins) move and sample the
membrane-plane dimension. When the distance between a lipid and a particular protein residue crosses a predefined threshold
we call it a contact. Our goal, thus, is to measure these contacts, keep track of them, and summarize them (using different parameters or metrics).

The main challenge here lies in distinguishing highly interacting lipids from other, transiently-interacting, ones. `prolintpy` provides
an intuitive API to calculate such contact-based metrics.


### Features
When calculating contact-based metrics, `prolintpy` supports the following features:

<ol>
<li>Automate the calculation of contact-based metrics.</li>
<li>Fully customizable.</li>
<li>Convert contacts into pandas DataFrames.</li>
<li>Extend contact definitions.</li>
<li>Easy to save and export contacts.</li>
</ol>

The exact value of the threshold is irrelevant and `prolintpy` gives you the option to provide whatever value.

### Contact Calculation
The contact interface of `prolintpy` is very intuitive. It uses the lipid and protein topologies as inputs, with many options available
to customize. We first define a contacts object:

```python
# takes as input the proteins and lipids topologis
contacts = pl.ComputeContacts(t, proteins, lipids)
```

That's it! Now, `contacts` contains all the required information, which allows us to easily calculate contacts, and customize them
without much effort and without having to change the code:

```python
# Calculate contact-based metrics between the provided protein and all residues
result = contacts.compute_neighbors(t)

# Calculate contact-based metrics between the provided protein and residue 45
c45 = contacts.compute_neighbors(t, [45])
print (c45)
> {'Protein0': {0: {45: <prolintpy.LPContacts for residue 45>}}}

# Calculate contact-based metrics between the provided protein and the following residues:
# the first 10, residue 50, and residues in the range 320-380
residues = [*range(1, 11), 50, *range(320, 381)]  # residues numbering has to match topology (i.e. no residue 0)
result = contacts.compute_neighbors(t, residues)
```

By default, `prolintpy` uses a cutoff of 0.7 nm to determine contacts, but this can be changed by simply providing the `cutoff` argument
to the `compute_neighbors` function. You can also use the `group` argument to group lipids according to their headgroup type (e.g. POPC,
PIPC, **PC lipids would be grouped into PC lipids, whereas POPE, PIPE, **PE lipids would be grouped into PE lipids, etc.). This is useful since it
allows you to measure average properties of different lipids:

```python
# Contacts using a 0.5 nm contact distance
c = contacts.compute_neighbors(t, cutoff=0.5)

# Contacts using a grouping of lipids
c = contacts.compute_neighbors(t, grouped=True)

# contacts using only the ROH beads of cholesterol lipids
c = contacts.compute_neighbors(t, atom_names=['ROH'])
```

### `prolintpy` helper functions

The output of a `compute_neighbors()` call is a python dictionary:

```python
c = contacts.compute_neighbors(t, [45])
print (c)
> {'Protein0': {0: {45: <prolintpy.LPContacts for residue 45>}}}
```
The innermost dictionary shows a contact object for the input residue (here 45). There is a similar object for each protein copy in the system (here 0).
All copies are grouped according to their protein name (here Protein0). Getting the contacts out of the dictionary boils down to providing the appropriate
keys:

```python
# get the contacts
print (c['Protein0'][0][45].contacts)
> {'CHOL': [600000.0, 900000.0, 1200000.0, 300000.0, 600000.0, 300000.0, 300000.0, 600000.0, 600000.0, 300000.0, 300000.0]}
```

ProLint also implements several helper functions to aid in getting these contacts. For instance:

```python
# Get contacts for residue 45
pl.retrieve_contacts(c, 45)
pl.retrieve_contacts_flat(c, 45)
```


### Contacts DataFrame

One of the helper functions provided by `prolintpy` is `contacts_dataframe` which builds a pandas DataFrame for all contacts. This is useful,
since many of the visualization applications rely on this dataframe structure. Using this function is straightforward once you have the output of the
compute_neighbors call.

```python
result = contacts.compute_neighbors(t)
df = pl.contacts_dataframe(result, proteins, t)
print (df.head())
>
```

If you are calculating contacts at multiple cutoff distances, then you can use the following code to build the DataFrame structure which, too, is
recognized by `prolintpy` visualization apps:

```python
radii = [0.45, 0.55, 0.65, 0.75, 0.85]
for radius in radii:
    result = contacts.compute_neighbors(t, cutoff=radius)

    if radius == radii[0]:
        df = pl.contacts_dataframe(result, proteins, t, radius, resolution)
    else:
        df = df.append(pl.contacts_dataframe(result, proteins, t, radius, resolution))
```

### Definitions
ProLint will calculate several contact metrics (e.g. average contact duration, occupancy, etc.). The following
are the metrics supported by default:

| Key         | Name                    | Description          |
| ----------- | ------------------- --- | -------------------- |
| `mean`      | Mean_Duration           | The average duration of all contacts|
| `max`       | Longest_Duration        | The longest duration (averaged if more than 1 protein).|
| `sum`       | Sum_of_all_Contacts     | The total sum of all contacts.|
| `lnr`       | Lipid_Number            | The average number of lipids in contact (the total number of contacts normalized with respect to the number of frames).|
| `nlnr`      | Normalized_Lipid_Number | Lipid_Number normalized with respect to the number of different lipids (e.g. number of different cholesterols).|
| `occ`       | Occupancy               | For each frame, we give a value of 0 if no lipid of interest is in contact with the residue and 1 otherwise. <br>Occupancy is then: sum_of_all_1s / nr_of_frames.|
| `rt`        | Residence_Time          | Residence Time (e.g. see <a href="https://www.pnas.org/content/117/14/7803">here</a> and <a href="https://pubs.acs.org/doi/abs/10.1021/ja310577u">here</a>. )|


You can easily define your own custom metrics.
We are planning on expanding the functionality of `prolintpy` with regards to contact-based metrics. So, your input/feedback is appreciated.
