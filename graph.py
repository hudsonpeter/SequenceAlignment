import matplotlib.pyplot as plt

ps = [16, 64, 128, 256, 384, 512, 768, 1024, 1280, 1536, 2048, 2560, 3072, 3584, 3968]

# for basic
memory = [11764, 11820, 11904, 12348, 13380, 14232, 16388, 19688, 25640, 31224, 49332, 64300, 96476, 99000, 139532]
time =[0.084, 0.322, 1.054, 4.17, 8.265, 14.48, 31.80, 58.34, 102.22, 137.585, 278.59, 409.25, 618.854, 680.65, 983.09]

print(len(ps))
print(len(memory))
print(len(time))
# for dc
memory_dc = [11864, 11900, 11680, 11860, 11948, 11836, 11744, 12300, 12196, 12344, 12304, 12176, 12512, 12628, 12520]
time_dc =[0.077, 0.616, 2.13, 7.608, 17.07, 29.687, 68.419, 119.97, 190.08, 263.58, 486.998, 756.727, 1053.52, 1467.59, 1828.96]
print(len(memory_dc))
print(len(time_dc))
plt.title('Memory(KB) vs Problem Size (M + N)')

plt.plot(ps, memory, label = "Basic", marker='o', markerfacecolor='blue', markersize=5)
plt.plot(ps, memory_dc, label = "Efficient", marker='o', markerfacecolor='orange', markersize=5)

plt.xlabel("Problem Size (M + N)")
plt.ylabel("Memory (KB)")
plt.legend()

plt.savefig('Memory_ProblemSize.png')
plt.show()

plt.title('Time(ms) vs Problem Size (M + N)')

plt.plot(ps, time, label = "Basic", marker='o', markerfacecolor='blue', markersize=5)
plt.plot(ps, time_dc, label = "Efficient", marker='o', markerfacecolor='orange', markersize=5)

plt.xlabel("Problem Size (M + N)")
plt.ylabel("Time(ms)")
plt.legend()

plt.savefig('Time_ProblemSize.png')
plt.show()