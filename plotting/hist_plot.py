data = """-1|5
-11
-13
-14
-14
-15
-15
-15
-15
-15|3|5
-16
-16|6
-2
-22
-22
-23|27
-23|39
-24
-24
-24
-25
-25
-25|27
-29
-29
-3
-3
-3
-3
-3
-3
-4
-4
-4
-6
-6
-6
-6
-7
-9
-9
-9
-9|-11
-9|-6
0
1
1|6
10
12
15
16
17
18
19
19|5
2
2
2
2
22
24
24
24|6
25
26
26
27|13
28
29|6
3
30
32
33
33
34
35
36
36
39
39
4
40
41|-15
41|6
42
42|20
45
47
48
49
49|5
5
5
5
5
5
5
5
5|-27
5|16
5|43
6
6
6
6
6
6
6
6
6
6
6
6
6
6
6
6|-1
6|-6
6|25|29
6|29
6|39
7|20
8
8
9
9"""
import matplotlib.pyplot as plt

data_list = data.split("\n")
curated_list = []
for i in data_list:
    if "|" not in i:
        curated_list.append(int(i))
    else:
        x = i.split("|")
        for j in x:
            curated_list.append(int(j))

curated_list_set = list(set(curated_list))
data_dict = {}
for i in curated_list_set:
    data_dict[i] = curated_list.count(i)
    if i == 5 or i == 6:
        print(data_dict[i])
max_freq = max(data_dict.values())

fig = plt.figure(figsize=(16, 9))
plt.bar(data_dict.keys(), data_dict.values())
#plt.hist(curated_list, bins=len(curated_list_set))
plt.title("Distances between binding sites and their downstream TSSs distribution")
plt.xlabel(f"TSS distances")
plt.ylabel("Frequency")
plt.yticks(range(0, max_freq + 1))
plt.xticks(range(min(curated_list), max(curated_list) + 1, 3))
plt.grid(True)
fig.savefig(f"TSS_distances_distribution.png")