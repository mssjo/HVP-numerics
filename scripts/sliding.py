TAB = '\t'

def convert_sliding(filename, stride):

    with open(f'{filename}.dat', 'r') as infile, open(f'{filename}_sliding{stride}.dat', 'w') as outfile:
        def output_average():
            if index < stride//2:
                # print(f"In head ({index}), skipping")
                return
            elif index < stride:
                # print(f"In head ({index}), at {buffers[0][index-stride//2]}, buffer {buffers[1]}, outputting {sum(buffers[1][:index+1])/(index+1)}")
                print(f"""{
                    buffers[0][index-stride//2]}{TAB}{
                    TAB.join(str(sum(buf[:index+1])/(index+1)) for buf in buffers[1:])
                    }""", file=outfile)
            else:
                # print(f"In body ({index}), at {buffers[0][(index-stride//2)%stride]}, buffer {buffers[1]}, outputting {sum(buffers[1])/stride}")
                print(f"""{
                    buffers[0][(index-stride//2)%stride]}{TAB}{
                    TAB.join(str(sum(buf)/stride) for buf in buffers[1:])
                    }""", file=outfile)

        def drain_buffers():
            nonlocal index
            while index % stride:
                # print(f"In tail ({index}), at {buffers[0][(index-stride//2)%stride]}, buffer {buffers[1]}, outputting {sum(buffers[1][index%stride:])/(stride - index%stride)}")
                print(f"""{
                    buffers[0][(index-stride//2)%stride]}{TAB}{
                    TAB.join(str(sum(buf[index%stride:])/(stride - index%stride)) for buf in buffers[1:])
                    }""", file=outfile)
                index += 1
            index = 0

        buffers = []
        for n, line in enumerate(infile):
            # print(f"{n}: {line.strip()}")
            if n == 0 or not line.strip():
                if n == 0:
                    for _ in line.strip().split():
                        buffers.append([0]*stride)
                    index = 0
                else:
                    drain_buffers()

                print(line, file=outfile)
            else:
                for i, val in enumerate(line.strip().split()):
                    buffers[i][index%stride] = float(val)

                output_average()
                index += 1

        drain_buffers()

        print(f"Read input from {infile.name}")
        print(f"Wrote stride-{10} sliding average to {outfile.name}")

if __name__ == '__main__':
    import sys

    convert_sliding(
        sys.argv[1][:-4] if sys.argv[1].endswith('.dat') else sys.argv[1],
        int(sys.argv[2]))
