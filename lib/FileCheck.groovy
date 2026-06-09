class FileCheck {
    static boolean checkoutfile(filePath) {
        def file = new File(filePath.toString())
        return file.exists() && file.isFile() && file.size() > 0
    }
}
