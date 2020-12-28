function doPost(to, p) { // to:提交动作（action）,p:参数 var myForm =
    var myForm = document.createElement("form");
    myForm.method = "post";
    myForm.action = to;
    for (var i in p) {
        var myInput = document.createElement("input");
        myInput.setAttribute("name", i);
        myInput.setAttribute("value", p[i]);
        myForm.appendChild(myInput);
    }
    document.body.appendChild(myForm);
    myForm.submit();
    console.log('done');
    document.body.removeChild(myForm); // 提交后移除创建的form }
}