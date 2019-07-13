import os


do = os.system

# do("python mm_run_with_pulling_start.py setup/6e67B/6e67B --to first_6e67B/4_no_contact_2 -s 1e7 --platform CUDA --reportFrequency 4000 -f forces_setup_6e67B.py --subMode -1")
# do("python mm_run_with_pulling_start.py setup/6e67B/6e67B --to first_6e67B/5_with_beta -s 1e7 --platform CUDA --reportFrequency 4000 -f forces_setup_6e67B.py")
# pdb_list = ["5xpd", "3kp9", "4a2n", "5d91"]
# pdb_list = ["4a2n", "5d91"]
pdb_list = ["2xov_complete", "6e67A", "5xpd", "3kp9", "4a2n", "5d91", "2jo1"]
pdb_list = ["2xov_complete", "6e67A", "5xpd", "3kp9", "4a2n", "5d91"]
# pdb_list = ["3kp9"]
# pdb_list = ["2xov_complete", "6e67A", "5xpd", "4a2n", "5d91"]
# pdb_list = ["5d91", "2xov_complete", "3kp9", "4a2n", "5xpd"]
# for pdb in pdb_list:
#     do(f"python mm_run_with_pulling_start.py setup/{pdb}/{pdb} --to first_{pdb}/1_with_contact -s 1e7 --platform CUDA --reportFrequency 4000 -f forces_setup_{pdb}.py --subMode 3")

# for i in range(3):
#     # batch prediction
#     for pdb in pdb_list:
#         do(f"python mm_run_with_pulling_start.py setup/{pdb}/{pdb} --to second_small_batch/{pdb}/{i} -s 1e6 --platform CUDA --reportFrequency 4000 -f forces_setup_{pdb}.py --subMode 3")


# for i in range(3):
#     # batch prediction
#     for pdb in pdb_list:
#         do(f"python mm_run_with_pulling_start.py setup/{pdb}/{pdb} --to second_small_batch/{pdb}/{i}_no_contact -s 1e6 --platform CUDA --reportFrequency 4000 -f forces_setup_{pdb}.py --subMode -1")


# for i in range(3):
#     # batch prediction
#     for pdb in pdb_list:
#         do(f"python mm_run_with_pulling_start.py setup/{pdb}/{pdb} --to second_small_batch/{pdb}/{i}_longer -s 4e6 --platform CUDA --reportFrequency 4000 -f energy_forces/forces_setup_{pdb}.py --subMode 3")

# for i in range(3):
#     # batch prediction
#     for pdb in pdb_list:
#         do(f"python mm_run_with_pulling_start.py setup/{pdb}/{pdb} --to second_small_batch/{pdb}/{i}_middle -s 2e6 --tempEnd 300 --platform CUDA --reportFrequency 4000 -f energy_forces/forces_setup_{pdb}.py --subMode 3")


# for i in range(1, 6):
#     # batch prediction
#     for pdb in pdb_list:
#         do(f"python mm_run_with_pulling_start.py setup/{pdb}/{pdb} --to thrid/{pdb}/{i} -s 1e6 --tempEnd 300 --platform CUDA --reportFrequency 4000 -f energy_forces/forces_setup_{pdb}.py --subMode 3")


# for i in range(1):
#     # batch prediction
#     for pdb in pdb_list:
#         do(f"python mm_run.py setup/{pdb}/{pdb} --to native/{pdb} -s 1e6 --tempStart 300 --tempEnd 300 --platform CUDA --reportFrequency 4000 -f energy_forces/forces_setup_{pdb}.py --subMode 3")

# for i in range(5):
#     # batch prediction
#     for pdb in pdb_list:
#         do(f"python mm_run_with_pulling_start.py setup/{pdb}/{pdb} --to fourth/{pdb}/{i} -s 1e6 --tempEnd 300 --platform CUDA --reportFrequency 4000 -f energy_forces/forces_setup_{pdb}.py --subMode 3")

# for i in range(5):
#     # batch prediction
#     for pdb in pdb_list:
#         do(f"python mm_run_with_pulling_start.py setup/{pdb}/{pdb} --to fifth_with_er/{pdb}/{i} -s 1e6 --tempEnd 300 --platform CUDA --reportFrequency 4000 -f energy_forces_with_er/forces_setup_{pdb}.py --subMode 3")

# for i in range(5):
#     # batch prediction
#     for pdb in pdb_list:
#         do(f"python mm_run_with_pulling_start.py setup/{pdb}/{pdb} --to sixth_with_er/{pdb}/{i} -s 2e6 --tempStart 1000 --tempEnd 300 --platform CUDA --reportFrequency 4000 -f energy_forces_with_er/forces_setup_{pdb}.py --subMode 3")

# for i in range(1):
#     # batch prediction
#     for pdb in pdb_list:
#         do(f"python mm_run.py setup/{pdb}/{pdb} --to native_short/{pdb}/{i} -s 1e5 --tempStart 300 --tempEnd 300 --platform CUDA --reportFrequency 4000 -f energy_forces_with_er/forces_setup_{pdb}.py --subMode 3")


# for i in range(20):
#     # batch prediction
#     for pdb in pdb_list:
#         do(f"python mm_run_with_pulling_start.py setup/{pdb}/{pdb} --to first/{pdb}/{i} -s 1e7 --tempStart 1000 --tempEnd 300 --platform CUDA --reportFrequency 4000 -f energy_forces_with_er/forces_setup_{pdb}.py --subMode 3")
#         do(f"python mm_run_with_pulling_start.py setup/{pdb}/{pdb} --to second_withoutExclusion/{pdb}/{i} -s 1e7 --tempStart 1000 --tempEnd 300 --platform CUDA --reportFrequency 4000 -f energy_forces_with_er_withoutExclusion/forces_setup_{pdb}.py --subMode 3")


for i in range(20):
    # batch prediction
    for pdb in pdb_list:
        do(f"python mm_run_with_pulling_start.py setup/{pdb}/{pdb} --to third_without_contact/{pdb}/{i} -s 1e7 --tempStart 1000 --tempEnd 300 --platform CUDA --reportFrequency 4000 -f energy_forces_with_er/forces_setup_{pdb}.py --subMode -10")



