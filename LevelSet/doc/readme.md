# �g����
## ���͉摜�̕��ƍ�����m��

�������s����B
```
./LevelSet.exe \
       --input /c/data/levelset/inputs/dreamworks.jpg \
       --mode info
```
�����̈Ӗ��͈ȉ��̒ʂ�
- input: ���͉摜�ւ̃p�X�ijpeg�̃O���C�摜�̂݃T�|�[�g�j
- mode: ��������m��Ƃ���info���w�肷��B

�����o�͂����B
```
  width :   820
 height :   376
```

## �����~��`�悷��
�������s����B
```
./LevelSet.exe \
        --input /c/data/levelset/inputs/dreamworks.jpg \
        --output /c/data/levelset/outputs \
        --mode initial \
        --x 494 \
        --y 164 \
        --radius 12 \
        --color 128
```
�����̈Ӗ��͈ȉ��̒ʂ�
- input: ���͉摜�ւ̃p�X�ijpeg�̃O���C�摜�̂݃T�|�[�g�j
- output: �o�̓t�H���_�ւ̃p�X
- mode: �����~��`�悷��Ƃ���initial���w�肷��B
- x: �~�̒��S��x���W
- y: �~�̒��S��y���W
- radius: �~�̔��a
- color: �~��`�悷��ۂ̐F�B[0,255]�͈̔͂Ō��߂�B
 
���s�̂���output�t�H���_�̒��Ɉȉ�2�̃t�@�C�����ł���B

- initial_loop.jpg: �����~���`�悳�ꂽ�摜
- circle_parameters.txt: �~�̒��S(x,y)�Ɣ��ar�ƐF���L�ڂ��ꂽ�e�L�X�g�t�@�C��

## LevelSet�@�����s����
�������s����B
```
./LevelSet.exe \
        --input /c/data/levelset/inputs/dreamworks.jpg \
        --output /c/data/levelset/outputs \
        --mode final \
        --x 494 \
        --y 164 \
        --radius 12 \
        --color 128 \
        --config /c/projects/levelset/x64/Release/levelset_config.txt \
        --interval 10
```
�����̈Ӗ��͈ȉ��̒ʂ�B
- input: ���͉摜�ւ̃p�X�ijpeg�̃O���C�摜�̂݃T�|�[�g�j
- output: �o�̓t�H���_�ւ̃p�X
- mode: �����~��`�悷��Ƃ���initial���w�肷��B
- x: �~�̒��S��x���W
- y: �~�̒��S��y���W
- radius: �~�̔��a
- color: �~��`�悷��ۂ̐F�B[0,255]�͈̔͂Ō��߂�B
- config: ��q����ݒ�t�@�C���ւ̃p�X
- interval: �摜��ۑ�����Ԋu

�ݒ�t�@�C�����̒��g�͈ȉ��̒ʂ�

![image](./config_sample.jpg)

�e�s�̃t�H�[�}�b�g��
```
<label><space><value>
```
�ł���B�e�p�����[�^�̈Ӗ��͈ȉ��̒ʂ�B

- time_step: ���Ԕ��W�̍��ݕ�
- time_step_number: ���Ԕ��W�̌J��Ԃ���
- space_step: ��ԍ��ݕ��Bx�������Ay�������Ƃ��ɓ����B
- constant_speed: ���x�֐��i��q�j�̒萔�����̒l�B�~������Ɏ��k������Ȃ畉�̒l���A���̋t�Ȃ琳�̒l��ݒ肷��B
- epsilon: ���x�֐��Ɍ����Â̒l�i��q�j
- sigma: ���͉摜�Ɏ{���K�E�X�ڂ����̓x�����B

���s�����output�t�H���_���Ɏ��X���X�ω�����Ȑ���`�悵���摜���o�͂����B��̏ꍇ�A���Ԕ��W��1000�X�e�b�v�itime_step_number�j�s���A�X�e�b�v10�iinterval�j���Ƃɉ摜���ۑ������B

## �⑫�F���x�֐��ɂ���

���x�֐��͎����Œ�`�����B

![image](./speed_function.jpg)

��1����contant_speed�ɁA��2���̃Â�epsilon�ɑ�������BK�͋ȗ���\���B
