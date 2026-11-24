# cpaudit_v1 | Colab Pro (A100) | upload cpaudit_v1.zip AND cpaudit_run_v1.py to
# MyDrive/NOMAD_Uncer_Quan/code_bundles/   then use one cell per step:

# Cell 1  (once: self-test, Phase 0, pre-registration; ~5-10 min)
from google.colab import drive; drive.mount("/content/drive")
%run "/content/drive/MyDrive/NOMAD_Uncer_Quan/code_bundles/cpaudit_run_v1.py" setup

# Cell 2  (Phase 1; long; re-run the same cell after a disconnect and it resumes)
from google.colab import drive; drive.mount("/content/drive")
%run "/content/drive/MyDrive/NOMAD_Uncer_Quan/code_bundles/cpaudit_run_v1.py" phase1

#   or in chunks:  ... cpaudit_run_v1.py" phase1 0-4     then 5-9, 10-14, 15-19, 20-24

# Cell 3  (report)
from google.colab import drive; drive.mount("/content/drive")
%run "/content/drive/MyDrive/NOMAD_Uncer_Quan/code_bundles/cpaudit_run_v1.py" report
