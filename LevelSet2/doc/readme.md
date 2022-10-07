# �g����
## �����Ȃ��Ŏ��s����

�������s����B
```
./LevelSet2.exe
```
�����o�͂����B
```
allowed options:
  --help                produce help message
  --dim arg             set either 2 or 3
  --verbose             print verbose description
  --input arg           set GRAY image path in 2D/set pattern in 3D
  --wband arg           set width of band
  --wreset arg          set width to reset
  --time_step arg       set time step
  --gain arg            set gain
  --constant_speed arg  set constant speed
  --speed_threshold arg set speed threshold
  --left arg            set left of initial rectangle
  --top arg             set top of initial rectangle
  --right arg           set right of initial rectangle
  --bottom arg          set bottom of initial rectangle
  --front arg           set front of initial rectangle (not used in 2D)
  --back arg            set back of initial rectangle (not used in 2D)
```
�����̈Ӗ��͈ȉ��̒ʂ�B
- dim: �Ώە��̎������w�肷��B�摜�Ȃ�2�A���̕��Ȃ�3�ł���B
- verbose: �f�o�b�O�����o�͂���B�l�������K�v�̂Ȃ������ł���B
- input: �摜�Ȃ�O���C�摜�ijpeg�̂݃T�|�[�g�j�ւ̃p�X���A���̕��Ȃ�STL�ւ̃p�X�������B
- wband: Narrow Band Level Set�@���̗p���Ă���B�g�ʁi���E�ʁj�̌v�Z�̈�̕���ݒ肷��i���}�Q�Ɓj�B
- wreset: Level Set�@�̐��x�����߂邽�߁A����v�Z�X�e�b�v���Ƃɔg�ʂ̏��������s���܂��B���̈������ƂȂ�̈�̕���ݒ肷��i���}�Q�Ɓj�B
- time_step: Level Set�@�ł͔g�ʂ����Ԕ��W������B���̎��ԕ���ݒ肷��B
- gain: �g�ʂ̑��x�͒萔���Ɣg�ʂ̋ȗ��Ɉˑ����鍀����Ȃ�B��҂̊�����ݒ肷��i���}�Q�Ɓj�B
- constant_speed: �g�ʑ��x�̒萔���i���}�Q�Ɓj
- speed_threshold: ���̑��x�ȉ��̏ꍇ�A�g�ʑ��x��0�Ƃ݂Ȃ��B
- left: �ŏ��ɗ^�����`���邢�͒����̂̃T�C�Y���w�肷��B
- top : �ŏ��ɗ^�����`���邢�͒����̂̃T�C�Y���w�肷��B
- right: �ŏ��ɗ^�����`���邢�͒����̂̃T�C�Y���w�肷��B
- bottom: �ŏ��ɗ^�����`���邢�͒����̂̃T�C�Y���w�肷��B
- front: �ŏ��ɗ^�����`���邢�͒����̂̃T�C�Y���w�肷��B3D�̏ꍇ�̂݁B
- back: �ŏ��ɗ^�����`���邢�͒����̂̃T�C�Y���w�肷��B3D�̏ꍇ�̂݁B
![image](./lsm.png)

��}�ɂ����āAwband�̂����΂Ŏ������̈�̕���wreset�ł���B�g�ʂ��΂̗̈�ɓ��B����ƁA�g�ʂ̍ď��������s�����B

���x�֐��͎����Œ�`�����B

![image](./speed_function.jpg)

��1����contant_speed�ɁA��2���̃Â�gain�ɑ�������B��2���͋Ȗʂ����炩�ɂ�����ʂ����BK�͋ȗ���\���B

## 3D�̏ꍇ
�������s����B
```
./LevelSet2.exe \
        --dim 3 \
        --verbose \
        --input "C:\data\arigis_datas\fr_vectors_3M\normalized_stls\S6617_DAK0319---.stl" \
        --wband 5 \
        --wreset 2 \
        --time_step 1  \
        --gain 0.1 \
        --constant_speed -1 \
        --speed_threshold 0.05 \
        --left 0 \
        --top 0 \
        --front 10 \
        --right 190 \
        --bottom 190 \
        --back 190
```
�E�B���h�E���N����OpenGL�ɂ��3D�`�󂪕`�悳���B�L�[�{�[�h���玟�̑��삪�\�ɂȂ�B
- ESC: �I�����܂��B�o�O���Ă�H�I���Ȃ��B
- p: �ꎞ��~���܂��B
- f: �g��(front)�̕`���L���E�����ɂ��܂��B
- o: ����(object)�̕`���L���E�����ɂ��܂��B
- <: ���_�𕨑̂ɋ߂Â��܂��B
- \>: ���_�𕨑̂��牓�����܂��B
- u: x������ɉ�]���܂��i�E�˂��j�B
- n: x������ɋt��]���܂��B
- h: y������ɉ�]���܂��i�E�˂��j�B
- j: y������ɋt��]���܂��B