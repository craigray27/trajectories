x=[0.2425782989
0.01346957451
0.38313885
0.4146526905
0.06776897286
0.9931269297
0.4843080465
0.765337766
0.0318338154
0.03093548167
0.9326404403
0.8878795341
0.5913297304
0.4787786847
0.8333543366
0.1863351968
0.7356527074
0.1150531718
0.6986586306
0.3556041114
0.6383000066
0.9082105076
0.2940004148
0.2649715889
0.3774939977
0.5416201155
0.009281815034
0.9994652765
0.01290191571
0.8424973208
0.8524702214
0.4670109839
0.05360723196
0.9767475906
0.1967547434
0.8569729407
0.1442137599
0.8006621882
0.7293976195
0.9857906471
0.1834057095
0.4997602107
0.4698620213
0.9709912995
0.4507702619
0.09579187729
0.9740816466
0.390235061
0.6806703558
0.0266693188
0.2312410745
0.4687390944
0.09796010055
0.4154100252
0.7962941871
0.3164026022
0.7785344495
0.8284933147
0.4871404206
0.36904878
0.6028454968
0.02426436079
0.8111117863
0.3557925929];
r=[0
3
6
9
12
15
18
21
23
26
29
32
35
38
41
44
46
49
52
55
58
61
64
67
69
72
75
78
81
84
87
90
92
95
98
101
104
107
110
113
115
118
121
124
127
130
133
136
138
141
144
147
150
153
156
159
161
163
165
167
169
171
173
175
176];
c=[0
1
8
1
2
9
2
3
10
3
4
11
4
5
12
5
6
13
6
7
14
7
15
8
9
16
9
10
17
10
11
18
11
12
19
12
13
20
13
14
21
14
15
22
15
23
16
17
24
17
18
25
18
19
26
19
20
27
20
21
28
21
22
29
22
23
30
23
31
24
25
32
25
26
33
26
27
34
27
28
35
28
29
36
29
30
37
30
31
38
31
39
32
33
40
33
34
41
34
35
42
35
36
43
36
37
44
37
38
45
38
39
46
39
47
40
41
48
41
42
49
42
43
50
43
44
51
44
45
52
45
46
53
46
47
54
47
55
48
49
56
49
50
57
50
51
58
51
52
59
52
53
60
53
54
61
54
55
62
55
63
56
57
57
58
58
59
59
60
60
61
61
62
62
63
63];
m=[325
-81
-81
325
-81
-81
325
-81
-81
325
-81
-81
325
-81
-81
325
-81
-81
325
-81
-81
325
-81
325
-81
-81
325
-81
-81
325
-81
-81
325
-81
-81
325
-81
-81
325
-81
-81
325
-81
-81
325
-81
325
-81
-81
325
-81
-81
325
-81
-81
325
-81
-81
325
-81
-81
325
-81
-81
325
-81
-81
325
-81
325
-81
-81
325
-81
-81
325
-81
-81
325
-81
-81
325
-81
-81
325
-81
-81
325
-81
-81
325
-81
325
-81
-81
325
-81
-81
325
-81
-81
325
-81
-81
325
-81
-81
325
-81
-81
325
-81
-81
325
-81
325
-81
-81
325
-81
-81
325
-81
-81
325
-81
-81
325
-81
-81
325
-81
-81
325
-81
-81
325
-81
325
-81
-81
325
-81
-81
325
-81
-81
325
-81
-81
325
-81
-81
325
-81
-81
325
-81
-81
325
-81
325
-81
325
-81
325
-81
325
-81
325
-81
325
-81
325
-81
325];
A=zeros(64,64);
c=c+1;
r=r+1;
q=1;
k=1;
while(k<=64)
for i=r(q):r(q+1)-1
    col=c(i);
    A(k,col)=m(i);
    if(col~=k)
        A(col,k)=m(i);
    end
end
    q=q+1;
    k=k+1;
end
p=zeros(64,64);
for i=1:64
    p(i,i)=1/A(i,i);
end
Z=p*A;
b=A*x;
xhat=A\b;
xhat2=[0.242578 0.0134696 0.383139 0.414653 0.067769 0.993127 0.484308 0.765338 0.0318338 0.0309355 0.93264 0.88788 0.59133 0.478779 0.833354 0.186335 0.735653 0.115053 0.698659 0.355604 0.6383 0.908211 0.294 0.264972 0.377494 0.54162 0.00928181 0.999465 0.0129019 0.842497 0.85247 0.467011 0.0536072 0.976748 0.196755 0.856973 0.144214 0.800662 0.729398 0.985791 0.183406 0.49976 0.469862 0.970991 0.45077 0.0957919 0.974082 0.390235 0.68067 0.0266693 0.231241 0.468739 0.0979601 0.41541 0.796294 0.316403 0.778534 0.828493 0.48714 0.369049 0.602845 0.0242644 0.811112 0.355793];
w=0;
matrix=zeros(1,288);
for i=1:64
    for j=1:64
        if(A(i,j)~=0)
            w=w+1;
            matrix(w)=A(i,j);
        end
    end
end

matrix2=[325 -81 -81 -81 325 -81 -81 -81 325 -81 -81 -81 325 -81 -81 -81 325 -81 -81 -81 325 -81 -81 -81 325 -81 -81 -81 325 -81 -81 325 -81 -81 -81 -81 325 -81 -81 -81 -81 325 -81 -81 -81 -81 325 -81 -81 -81 -81 325 -81 -81 -81 -81 325 -81 -81 -81 -81 325 -81 -81 -81 -81 325 -81 -81 325 -81 -81 -81 -81 325 -81 -81 -81 -81 325 -81 -81 -81 -81 325 -81 -81 -81 -81 325 -81 -81 -81 -81 325 -81 -81 -81 -81 325 -81 -81 -81 -81 325 -81 -81 325 -81 -81 -81 -81 325 -81 -81 -81 -81 325 -81 -81 -81 -81 325 -81 -81 -81 -81 325 -81 -81 -81 -81 325 -81 -81 -81 -81 325 -81 -81 -81 -81 325 -81 -81 325 -81 -81 -81 -81 325 -81 -81 -81 -81 325 -81 -81 -81 -81 325 -81 -81 -81 -81 325 -81 -81 -81 -81 325 -81 -81 -81 -81 325 -81 -81 -81 -81 325 -81 -81 325 -81 -81 -81 -81 325 -81 -81 -81 -81 325 -81 -81 -81 -81 325 -81 -81 -81 -81 325 -81 -81 -81 -81 325 -81 -81 -81 -81 325 -81 -81 -81 -81 325 -81 -81 325 -81 -81 -81 -81 325 -81 -81 -81 -81 325 -81 -81 -81 -81 325 -81 -81 -81 -81 325 -81 -81 -81 -81 325 -81 -81 -81 -81 325 -81 -81 -81 -81 325 -81 -81 325 -81 -81 -81 325 -81 -81 -81 325 -81 -81 -81 325 -81 -81 -81 325 -81 -81 -81 325 -81 -81 -81 325 -81 -81 325 -81];
matrix-matrix2;

