# fuzzy-match
## How the core similarity check works:
Represent the string as a set of bit arrays for each character in the alphabet, this is why this works best for small alphabets like DNA.
For example the string `AGCTCCGT` would be represented as four arrays:
```
A = 10000000
C = 00101100
G = 01000010
T = 00010001
```

These arrays are then stored in a contiguous memory so for example this string would be represented as the 32b number `10000000 00101100 01000010 00010001`.
This representations is called DNA4 

To check the number of matched characters between two strings we AND them together and then count the set bits.
For example the number of matches between the above string and ACGCTCGA could be calculated as so.
```
ACGCTCGA -> 10000001 01010100 00100010 00001000
A = 10000001
C = 01010100
G = 00100010
T = 00001000
```
```
AGCTCCGT & ACGCTCGA = 
10000000 00101100 01000010 00010001
&
10000001 01010100 00100010 00001000
=
10000001 00000100 00000010 00000000
```
We can then simply call a bit counting instruction like popcount to find the number of matches is 4.
Creating a DNA4
### Each character has a bit pattern associated with it
```
A = 00000001 00000000 00000000 00000000
C = 00000000 00000001 00000000 00000000
G = 00000000 00000000 00000001 00000000
T = 00000000 00000000 00000000 00000001
```
Simply move along the string adding the relevant bit pattern and then shifting the result left by 1.

## Implementation
Currently the length is hard coded as 32 characters so a string is represented by two 64b ints. The ints are padded with 0's to make this length so the actual representation of the above sequence would be:
```
uint64_t AC = 00000000 00000000 00000000 10000000 00000000 00000000 00000000 00101100
             |         A's Padding      |   A's  |         C's Padding      |   C's  |

uint64_t GT = 00000000 00000000 00000000 01000010 00000000 00000000 00000000 00010001
             |         G's Padding      |   G's  |         T's Padding      |   T's  |
```
For DNA4 creation the character bit patterns are also grouped by AC and GT as 64b uints.
			 
Putting everything together to find the similarity between two strings a and b that have been converted to this form we perform the following 4 instructions:
```
AC_matches = a.AC & b.AC
GT_matches = a.GT & b.GT
all_matches = AC_matches | GT_matches
similarity = popcount(all_matches)
```
### Performance of Simple Solution:
This lets us check the similarity between two strings in just 4 operations but the strings need to be in DNA4 form which takes length(string) time to convert. For this reason this method is bad at matching small numbers of strings together. However in the case where you are matching a string to all locations within a long sequence (such as target matching) the trade-off is worth it.

From here on we'll call the strings (converted to DNA4) that are being matched to a sequence **targets** and the sequence itself **seq**.

Matching several targets to a sequence.
1. First all strings to match are converted to their DNA4 form.
2. Then the first length(string) characters of the sequence are converted to DNA4 form.

3. All targets are matched to seq to check similarity.
4. seq is moved along one character by adding the character to the end, and bit shifting left.
We do not not need to zero off the left side of seq because the targets will have zeros there so the AND will always be zero at those locations.

Repeat steps 3 and 4 until the sequence is consumed. 

If we let sequence length be S, number of targets be T and length of targets be L:
The total number of required instructions is around ST \* (4 + if check on similarity + loop overhead). (Cost of the IF is usually hidden as a rare branch)
This assumes that ST >> S+LT such that the cost of converting to DNA4 representation is amortized.

## Optimizations
Although this is fast it is still a brute force solution of checking every possible target against every possible position in the sequence. This is impractical with large datasets where S\*T may be orders of magnitude over a trillion.

Imagine a case where L = 30
In a perfect matching case we could simply store the targets hashed by say their first 10 characters then for each position in seq we only need to check the targets with the same leading hash. This would decrease the number of comparisons by a factor 4^10, ~ 1 millionth of the comparisons.
 
Consider instead the case where we allow a mismatch, so we may have a mismatch in the leading 10 characters. To allow for this we store the targets in 2 hashmaps, hashed by the first and second 10 characters. Any target that differs in no more than 1 location from the seq will have an exact match in at least 1 of the sections. This doubles the number of hash buckets we need to check but each bucket still contains a millionth of the number of targets, definitely a worthy trade-off.

However what if we allowed several mismatches, say, up to M where M < L.
We require M+1 segments to guarantee that one is a perfect match. The buckets must not overlap so the size of the buckets must be floor(L/(M+1)), for example if L is 30 and M is 7 the bucket size would be 3. In this case each bucket would have the 4^3 times fewer targets in them but we would need to check 10 buckets for a total speed up over the brute force solution of just over 6. 
So the total number of comparisons required would be S\*floor(L/(M+1))/(4^floor(L/(M+1))).

What if we allow mismatches within each division. If in the above example we instead divide the target into 4 segments then there necessarily must be a segment with 1 or fewer mismatches present. If we make 7 hash maps for each segment, each with a different character left out then one of THOSE hashes must be a perfect match, the hash consists of 6 characters so only contains T/(4^6) targets. This leads to checking 28 total hashmaps but each bucket only contains T/4096 targets. A speed up of 146 over the brute force solution.

There is nothing special about allowing only 1 mismatch within each division. In the general case if we split L into N divisions (N<=L) then the size of each division D is floor(L/N), the minimum number of mismatches present in a division K is floor(M/N), and so the target reduction factor is 4^(D-K), the total number of hashmaps is the number of divisions times the number of ways the required mismatches can be removed (D choose K).
So the total average number of targets to check is calculated as following:

T_compared = T\*N\*(floor(L/N) choose floor(M/N))/(4^(floor(L/N)-floor(M/N)))

For a given run L and M are constant so to calculate the optimal value of N we simply test the above equation for all values of N<=L to find the minimum.

There are a couple of caveats to this however. Creating these hash maps isn't free, they take time and space. checking 100000 hash_maps with each bucket only containing a couple of targets may be fewer comparisons than 1 hash_map with a million targets in each bucket. However, the cost of all the redirection and cache misses will quickly eliminate this benefit. In addition the number of hashes that need to be calculated is (T+S)\*N\*(floor(L/N) choose floor(M/N)), and these hashes are far more expensive than comparisons (on the order of (floor(L/N)-floor(M/N))). For this reason the hash cost is also added to the calculation of the optimum value of N.
```
total runtime = O(S*T*N*(floor(L/N) choose floor(M/N))/(4^(floor(L/N)-floor(M/N)))+(T+S)*N*(floor(L/N) choose floor(M/N)) * (floor(L/N)-floor(M/N)))

total runtime = O(N*(floor(L/N) choose floor(M/N))((S*T)/(4^(floor(L/N)-floor(M/N)))+(T+S)*(floor(L/N)-floor(M/N))))
                 |        Number of hash maps    ||             Comparisons        | |           Hashing         |
```
Finally there is a space constraint, the total space of all these hash maps is:

```
Hash_space_in_Bytes = N*(floor(L/N) choose floor(M/N))*(16*T + 24*4^(floor(L/N)-floor(M/N))))
                     |        Number of hash maps    ||Targets|    Vectors for hashmaps     |
```
This can be extremely large if T is large. Consider the example where T = 1million, L = 30 and M is 7 but say we select N = 3 instead of N = 2:
There appears to be less than 30% of the comparisons as when N = 4 however the space requirement goes from 450 MB to 2.4 GB
For this reason we provide a size hint option for the maximum size to allow for the index.
