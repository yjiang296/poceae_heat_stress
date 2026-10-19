# Arial字体

请将自己有权使用的Arial TrueType字体复制为本目录的 `Arial.ttf`。仓库不分发字体二进制文件。

`fonts.conf` 使用相对于自身位置的字体目录，支持移动整个流程目录。批量提交器会将本机字体复制到批次代码快照，并生成快照专用字体配置，以保持每次作业的字体一致。

```bash
cp /path/to/your/Arial.ttf assets/Arial.ttf
FONTCONFIG_FILE="$PWD/assets/fonts.conf" fc-match Arial
```

在 `DEG_stringtie/` 目录运行。应返回Arial字体；没有Arial时先完成合法安装或复制，再运行绘图和提交器。步骤01—03本身不使用字体。
